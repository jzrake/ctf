
local optparse = require 'optparse'
local MaraSim  = require 'MaraSim'
local util     = require 'util'
local oo       = require 'class'
local sim      = require 'simulation'
local array    = require 'array'
local hdf5     = require 'lua-hdf5.LuaHDF5'
local unigrid  = require 'unigrid'
local cow      = require 'cow'
local MPI      = require 'MPI'
local Mara     = require 'Mara'
local problems = require 'problems'

local MyBase = oo.class('MyBase', sim.SimulationBase)
local MyMara = oo.class('MyMara', MaraSim.MaraSimulation)


function MyMara:initialize_physics()
   local function pinit(x,y,z)
      return self.problem:solution(x,y,z,0)
   end
   local P = self.primitive.array
   Mara.init_prim(P:buffer(), pinit)
end

function MyMara:initialize_solver()
   local opts = self.user_opts

   self.CFL = opts.CFL or 0.4
   self.Ng = 3
   self.Nx = opts.resolution or 128
   self.Ny = opts.resolution or 128
   self.Nz = 1
   self.dx = 1.0 / self.Nx
   self.dy = 1.0 / self.Ny


   MPI.Init()
   cow.init(0, nil, 0) -- to reopen stdout to dev/null

   local domain = cow.domain_new()
   local domain_comm = MPI.Comm()

   cow.domain_setndim(domain, 2)
   cow.domain_setsize(domain, 0, self.Nx)
   cow.domain_setsize(domain, 1, self.Ny)
   cow.domain_setsize(domain, 2, self.Nz)
   cow.domain_setguard(domain, self.Ng)
   cow.domain_commit(domain)
   cow.domain_getcomm(domain, domain_comm)

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
   Mara.set_advance(advance or 'rk3')
   Mara.set_godunov(solver)
   Mara.set_boundary(self.problem:boundary_conditions())
   Mara.set_riemann(self.user_opts.riemann or 'hllc')
   Mara.config_solver({theta=2.0}, true)

   local prim_names = Mara.fluid.GetPrimNames()
   local Nq = #prim_names

   self.domain = domain
   self.primitive = unigrid.DataManagerHDF5(domain, prim_names)
   self.shape = self.primitive:local_mesh_size('shape')
   Mara.set_domain({0,0}, {1,1}, self.shape, Nq, self.Ng, domain_comm)
end

function MyMara:initialize_behavior()
   local opts = self.user_opts
   local cpi = opts.cpi or 0.1
   local tmax = opts.tmax or self.problem:finish_time()
   self.behavior.message_cadence = opts.message_cadence or 1
   self.behavior.checkpoint_cadence = tonumber(cpi)
   self.behavior.max_simulation_time = tonumber(tmax)
end

function MyMara:advance_physics()
   local dt = self.status.time_increment
   local P = self.primitive.array:buffer()
   Mara.advance(P, dt)
end

function MyMara:local_mesh_size()
   return self.primitive:local_mesh_size()
end

function MyMara:checkpoint_write(fname)
   local base = self.user_opts.id or 'chkpt'
   local n = self.status.checkpoint_number
   local fname = fname or string.format('data/%s.%04d.h5', base, n)
   self.primitive:write(fname, {group='prim'})
end

function MyMara:user_work_finish()
   self.problem:user_work_finish()
end

function MyMara:user_work_iteration()
   self.problem:user_work_iteration()
end

function MyMara:finalize_solver()
   cow.domain_del(self.domain)
   Mara.close()
   MPI.Finalize()
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
   parser.add_option{"--plot", dest="plot", action="store_true"}
   parser.add_option{"--reconstruction", dest="reconstruction"}
   parser.add_option{"--riemann", dest="riemann"}
   parser.add_option{"--advance", dest="advance",
                     help="which Runge-Kutta to use for solution advance"}
   parser.add_option{"--solver", dest="solver"}
   parser.add_option{"--resolution", "-N", dest="resolution", help="grid resolution"}
   parser.add_option{"--message-cadence", dest="message_cadence",
                     help="print a message every N iterations"}

   local opts, args = parser.parse_args()
   local problem_class = problems[arg[2]]

   local sim = MyMara(opts)
   local problem = problems.SmoothKelvinHelmholtz(opts)
   sim:run(problem)
end

main()
