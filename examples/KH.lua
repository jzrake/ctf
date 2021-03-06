
local optparse = require 'optparse'
local MaraSim  = require 'MaraSim'
local util     = require 'util'
local oo       = require 'class'
local sim      = require 'simulation'
local array    = require 'array'
local hdf5     = require 'lua-hdf5.LuaHDF5'
local unigrid  = require 'new-unigrid'
local cow      = require 'cow'
local MPI      = require 'MPI'
local Mara     = require 'Mara'
local problems = require 'problems'
local json     = require 'json'

local MyBase = oo.class('MyBase', sim.SimulationBase)
local MyMara = oo.class('MyMara', MaraSim.MaraSimulation)

function MyMara:initialize_physics()
   local function pinit(x,y,z)
      return self.problem:solution(x,y,z,0)
   end
   local P = self.Primitive
   Mara.init_prim(P:buffer(), pinit)
end
function MyMara:domain_dimensions() return self.ndim end
function MyMara:initialize_solver()
   local opts = self.user_opts

   self.CFL = opts.CFL or 0.4
   self.Ng = 3
   self.Nx = opts.resolution or 32
   self.Ny = opts.resolution or 32
   self.Nz = opts.resolution or 32
   self.ndim = tonumber(opts.ndim) or 3

   MPI.Init()
   cow.init(0, nil, 0) -- to reopen stdout to dev/null

   local fluid = ({nrhyd='euler',
                   srhyd='srhd',
                   srmhd='rmhd'})[self.problem:fluid()]
   if not fluid then
      error('Mara does not support fluid system '..self.problem:fluid())
   end

   local advance = self.user_opts.advance
   local solver = ({
                      spectral = 'weno-split',
                      godunov = 'plm-split',
                      muscl = 'plm-muscl'})[self.user_opts.solver or 'godunov']

   if self.user_opts.solver == 'muscl' and advance ~= 'single' then
      print('[MaraSim] Warning! --solver=muscl only supports --advance=single'
         ..', going with single')
         advance = 'single'
   end

   Mara.start()
   Mara.set_fluid(fluid)
   Mara.set_eos('gamma-law', self.problem:adiabatic_index())
   Mara.set_advance(advance or 'rk3')
   Mara.set_godunov(solver)
   Mara.set_boundary(self.problem:boundary_conditions())
   Mara.set_riemann(self.user_opts.riemann or 'hllc')
   Mara.config_solver({theta  =  opts.plm_theta or 2.0,
                       IS     =  opts.IS or 'js96',
                       sz10A  =  opts.sz10A or 100.0}, false)

   local prim_names = Mara.fluid.GetPrimNames()
   local Nq = #prim_names
   local Ng = self.Ng
   local L0, L1 = self.problem:domain_extent()
   local N = { }
   for d=1,3 do
      if self.ndim >= d then
         N[d] = ({self.Nx, self.Ny, self.Nz})[d]
      else
	 L0[d] = nil
	 L1[d] = nil
      end
   end

   local domain = unigrid.UnigridDomain(N, Ng)
   local prim_manager = unigrid.UnigridDataField(domain, prim_names)
   local domain_comm = domain:get_comm()

   self.cart_rank = domain:get_rank()
   self.Primitive = prim_manager:array()
   self.domain = domain
   self.prim_manager = prim_manager
   Mara.set_domain(L0, L1, N, Nq, Ng, domain_comm)

   if hdf5.have_mpio() then
      domain:set_collective(1)
   end
end

function MyMara:initialize_behavior()
   local opts = self.user_opts
   local cpi = opts.cpi or 0.1
   local tmax = opts.tmax or self.problem:finish_time()
   self.behavior.run_identifier = opts.id or 'test'
   self.behavior.message_cadence = opts.message_cadence or 1
   self.behavior.checkpoint_cadence = tonumber(cpi)
   self.behavior.max_simulation_time = tonumber(tmax)
   self.measure_log = { }
end

function MyMara:local_mesh_size()
   return self.domain:size('all', 'local-interior')
end

function MyMara:checkpoint_write(fname)
   local base = 'chkpt'
   local id = self.behavior.run_identifier
   local n = self.status.checkpoint_number
   local fname = fname or string.format('data/%s/%s.%04d.h5', id, base, n)

   if self.cart_rank == 0 then
      if util.file_exists(fname) then
	 os.execute("rm "..fname)
      end
      os.execute('mkdir -p data/'..id)
   end

   self.prim_manager:write(fname, 'prim')
   self.domain:barrier()
   if self.cart_rank == 0 then
      local chkpt = hdf5.File(fname, 'r+')
      chkpt["measure_log"] = json.encode(self.measure_log)
      chkpt["status"] = json.encode(self.status)
      if self.problem.model_parameters then
	 chkpt["model_parameters"] = json.encode(self.problem.model_parameters)
      end
      chkpt:close()
   end
end

function MyMara:user_work_finish()
   local output = self.user_opts.output
   if output then
      self:checkpoint_write(output)
   end
   self.problem:user_work_finish()
end


-- *****************************************************************************
-- Function to collect all available measurements from Mara
-- .............................................................................
function MyMara:user_work_iteration()
   self.problem:user_work_iteration()
   if self.status.iteration_number % 4 == 0 then
      local meas = { }
      meas.U                    = Mara.measure.mean_cons()
      meas.P                    = Mara.measure.mean_prim()
      meas.energies             = Mara.measure.mean_energies()
      meas.mean_velocity        = Mara.measure.mean_velocity()
      meas.mean_T, meas.max_T   = Mara.measure.mean_max_temperature()
      meas.mean_B, meas.max_B   = Mara.measure.mean_max_magnetic_field()
      meas.mean_Ms, meas.max_Ms = Mara.measure.mean_max_sonic_mach()
      meas.mean_Ma, meas.min_Ma = Mara.measure.mean_min_alfvenic_mach()
      meas.max_lorentz_factor   = Mara.measure.max_lorentz_factor()
      meas.status = util.deepcopy(self.status)
      table.insert(self.measure_log, meas)
   end
end

function MyMara:finalize_solver()
   Mara.close()
   self.domain = nil
   self.prim_manager = nil
   collectgarbage()
   MPI.Finalize()
end

function handle_crash_srmhd(self, attempt)
   local P = self.Primitive:buffer()
   local status = self.status
   local r = 0.0
   Mara.set_advance("single")
   if attempt == 0 then -- healthy time-step
      Mara.set_godunov("plm-muscl")
      Mara.set_riemann("hlld")
      Mara.config_solver({theta=2.0}, true)
      return 0
   elseif attempt == 1 then
      status.time_increment = 0.5 * status.time_increment
      return 0
   elseif attempt == 2 then
      status.time_increment = 0.5 * status.time_increment
      return 0
   elseif attempt == 3 then
      Mara.set_riemann("hll")
      status.time_increment = 0.5 * status.time_increment
      return 0
   elseif attempt == 4 then
      Mara.config_solver({theta=1.0}, true)
      status.time_increment = 0.5 * status.time_increment
      return 0
   elseif attempt == 5 then
      Mara.config_solver({theta=0.0}, true)
      status.time_increment = 0.5 * status.time_increment
      return 0
   else
      return 1
   end
end


local function main()
   local usage = "tests-1d <problem> [<options>]"
   local parser = optparse.OptionParser{usage=usage,
                                        version="CTF version 1.0"}

   parser.add_option{"--cpi", dest="cpi", help="checkpoint interval"}
   parser.add_option{"--id", dest="id", help="problem ID: used for checkpoint names"}
   parser.add_option{"--cfl", dest="CFL",
                     help="Courant-Freidrichs-Lewy time-step constraint"}
   parser.add_option{"--tmax", dest="tmax", help="end simulation time"}
   parser.add_option{"--plot", dest="plot", action="store_true",
                     help="launch a gnuplot window after run (1d only)"}
   parser.add_option{"--reconstruction", dest="reconstruction"}
   parser.add_option{"--IS", dest="IS",
                     help='WENO smoothness indicator: [js96, b08, sz10]'}
   parser.add_option{"--sz10A", dest="sz10A",
                     help="Shen & Zha (2010) 'A' parameter: [0,100]"}
   parser.add_option{"--reconstruction", dest="reconstruction"}
   parser.add_option{"--riemann", dest="riemann"}
   parser.add_option{"--advance", dest="advance",
                     help="which Runge-Kutta to use for solution advance"}
   parser.add_option{"--solver", dest="solver",
                     help="godunov (riemann solver) or spectral (characteristic-wise)"}
   parser.add_option{"--resolution", "-N", dest="resolution", help="grid resolution"}
   parser.add_option{"--ndim", "-d", dest="ndim", help="domain dimensionality"}
   parser.add_option{"--message-cadence", dest="message_cadence",
                     help="print a message every N iterations"}
   parser.add_option{"--output", "-o", dest="output",
                     help="write an HDF5 file of the final solution"}
   parser.add_option{"--model_parameters", "-p", dest="model_parameters",
		     help="extra problem-specific parameters (as a table)"}
   parser.add_option{"--problem", dest="problem", help="name of a problem class"}

   local opts, args = parser.parse_args()
   local problem_class = problems[opts.problem]

   if not problem_class then
      print("[!] provide a problem name with --problem=ProblemClass")
      print("[!] available problems are")
      util.pretty_print(problems)
      os.exit()
   end

   local sim = MyMara(opts)
   local problem = problem_class(opts)

   if oo.isinstance(problem, problems.TearingMode) then
      sim.handle_crash = handle_crash_srmhd
   end
   if oo.isinstance(problem, problems.MagneticBubble) then
      sim.handle_crash = handle_crash_srmhd
   end

   sim:run(problem)
end

main()
