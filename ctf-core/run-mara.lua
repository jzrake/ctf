
if arg[2] == '--explain' then
   print "run a 2d/3d parallel research problem with Mara"
   return
end

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
   --
   -- Custom boundary conditions configuration for MagneticBubble problem
   --
   local bc = self.problem:boundary_conditions()
   if bc == 'magnetic-bubble' then
      local r0 = self.problem.model_parameters.r0
      local prof = self.problem.model_parameters.prof
      Mara.set_boundary(bc, r0, prof)
   elseif oo.classname(self.problem) == 'Magnetar' then
      local L0 = self.problem.model_parameters.L0
      local B0 = self.problem.model_parameters.B0
      local C0 = self.problem.model_parameters.C0
      --Mara.set_fluxsrc('magnetar', L0, B0, C0)
      Mara.set_srcterm('magnetar', L0, B0, C0)
      Mara.config_solver({preset=self.problem.model_parameters.P0,
			  dreset=self.problem.model_parameters.D0}, false)
   elseif oo.classname(self.problem) == 'Wind' then
      Mara.set_srcterm('wind')
   end

   --
   -- Generic initial model data
   --
   local function pinit(x,y,z)
      return self.problem:solution(x,y,z,0)
   end
   local P = self.Primitive
   if self.user_opts.restart then
      self.prim_manager:read(self.user_opts.restart, 'prim')
      Mara.advance(P:buffer(), 0.0) -- copies array to Mara internal memory
   else
      Mara.init_prim(P:buffer(), pinit)
   end
end
function MyMara:domain_dimensions() return self.ndim end
function MyMara:initialize_solver()
   local opts = self.user_opts

   if opts.restart then
      local chkpt = hdf5.File(opts.restart, 'r')
      local extent = chkpt['prim']['rho']:get_space():get_extent()
      if opts.Nx or opts.Ny or opts.Nz or opts.resolution then
	 print('[run-mara] Warning! command-line option specified a domain'..
	       ' resolution which was ignored (inferred from checkpoint file)')
      end
      self.Nx = extent[1]
      self.Ny = extent[2]
      self.Nz = extent[3]
      self.ndim = #extent
   else
      self.Nx = opts.Nx or opts.resolution or 32
      self.Ny = opts.Ny or opts.resolution or 32
      self.Nz = opts.Nz or opts.resolution or 32
      self.ndim = tonumber(opts.ndim) or 3
   end

   self.Ng = 3
   self.CFL = opts.CFL or 0.4

   if not opts.serial then
      MPI.Init()
      cow.init(0, nil, 0) -- to reopen stdout to dev/null
   end

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
      print('[run-mara] Warning! --solver=muscl only supports --advance=single'
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
   print("domain extent:")
   print(unpack(L0))
   print(unpack(L1))
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

   if self.user_opts.restart then
      local chkpt = hdf5.File(opts.restart, 'r')
      local status_string = chkpt['status']:value()
      local measure_log_string = chkpt['measure_log']:value()
      self.status = json.decode(status_string)
      self.measure_log = json.decode(measure_log_string)
   end
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
      chkpt["git_sha"] = Mara.git_sha()
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
   if not self.user_opts.serial then
      MPI.Finalize()
   end
end

local handle_crash = { }
function handle_crash.Magnetar(self, attempt)
   local P = self.Primitive:buffer()
   local status = self.status
   local r = 0.5
   Mara.set_advance("rk3")
   if attempt == 0 then -- healthy time-step
      Mara.set_godunov("plm-split")
      Mara.set_riemann("hlld")
      Mara.config_solver({theta=2.0, pfloor=1e-6, ereset=false}, true)
      return 0
   elseif attempt == 1 then
      Mara.diffuse(P, 0.5)
      return 0
   elseif attempt == 2 then
      Mara.diffuse(P, 0.5)
      return 0
   elseif attempt == 3 then
      Mara.diffuse(P, 0.5)
      return 0
   elseif attempt == 4 then
      Mara.diffuse(P, 0.5)
      return 0
   elseif attempt == 5 then
      Mara.config_solver({theta=1.5}, true)
      return 0
   elseif attempt == 6 then
      Mara.config_solver({theta=1.0}, true)
      return 0
   elseif attempt == 7 then
      Mara.config_solver({theta=0.0, ereset=true}, true)
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
   parser.add_option{"--resolution", "-N", dest="resolution",
		     help="default grid resolution"}
   parser.add_option{"--Nx", dest="Nx", help="grid resolution along x-axis"}
   parser.add_option{"--Ny", dest="Ny", help="grid resolution along y-axis"}
   parser.add_option{"--Nz", dest="Nz", help="grid resolution along z-axis"}
   parser.add_option{"--ndim", "-d", dest="ndim", help="domain dimensionality"}
   parser.add_option{"--message-cadence", dest="message_cadence",
                     help="print a message every N iterations"}
   parser.add_option{"--output", "-o", dest="output",
                     help="write an HDF5 file of the final solution"}
   parser.add_option{"--model_parameters", "-p", dest="model_parameters",
		     help="extra problem-specific parameters (as a table)"}
   parser.add_option{"--problem", dest="problem", help="name of a problem class"}
   parser.add_option{"--restart", dest="restart",
		     help="restart from named checkpoint file"}
   parser.add_option{"--serial", "-s", dest="serial", action="store_true",
		     help="run a serial job, disable MPI initialization"}

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

   sim.handle_crash = handle_crash[oo.classname(problem)]
   sim:run(problem)
end

main()
