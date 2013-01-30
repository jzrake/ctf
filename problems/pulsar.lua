
local sim      = require 'simulation'
local oo       = require 'class'
local array    = require 'array'
local util     = require 'util'
local hdf5     = require 'lua-hdf5.LuaHDF5'
local Mara     = require 'Mara'
local problems = require 'problems'


local MaraSimulation = oo.class('MaraSimulation', sim.SimulationBase)
local NSAccretion = oo.class('NSAccretion', problems.TestProblem)

function MaraSimulation:initialize_behavior()
   local opts = self.user_opts
   local cpi = opts.cpi or 1.0
   local tmax = opts.tmax or 1.0
   local dynamical_time = self.problem:dynamical_time()
   self.behavior.message_cadence = 1
   self.behavior.checkpoint_cadence = cpi * dynamical_time
   self.behavior.max_simulation_time = tmax * dynamical_time
end

function MaraSimulation:initialize_solver()
   local opts = self.user_opts
   self.CFL = opts.CFL or 0.4
   self.Ng = 3
   self.N = opts.resolution or 128

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

   Mara.set_domain({-1, -1}, {1, 1}, {self.N, self.N}, Nq, self.Ng)
end

function MaraSimulation:report_configuration()
   Mara.show()
   Mara.units.Print()
end

function MaraSimulation:finalize_solver()
   Mara.close()
end

function MaraSimulation:initialize_physics()
   local P, G = self.problem:solution(0.0)
   self.Primitive = P
   self.Gravity = G
end

function MaraSimulation:set_time_increment()
   local dt = Mara.get_timestep(self.CFL)
   self.status.time_increment = dt == math.huge and 0.0 or dt
end

function MaraSimulation:advance_physics()
   local dt = self.status.time_increment
   local P = self.Primitive:buffer()
   local time, err = Mara.advance(P, dt)
   if err ~= 0 then
      error('Mara.advance returned error code '..err)
   end
end

function MaraSimulation:checkpoint_write()
   local n = self.status.checkpoint_number
   local t = self.status.simulation_time
   local Ng = self.Ng
   local fname = string.format('data/chkpt.%04d.h5', n)
   local outfile = hdf5.File(fname, 'w')
   outfile['prim' ] = self.Primitive[{{Ng,-Ng},{Ng,-Ng},nil}]
   outfile:close()
end

function MaraSimulation:user_work_iteration()
   self.problem:user_work_iteration()
end

function MaraSimulation:user_work_finish()
   local t = self.status.simulation_time
   local Ng = self.Ng
   local P = self.Primitive
   local Pexact = self.problem:solution(t)
   local dx = self.dx
   local P0 = Pexact[{{Ng,-Ng},nil}]:vector()
   local P1 = P     [{{Ng,-Ng},nil}]:vector()
   local L1 = 0.0
   for i=0,#P0/5-1,5 do
      L1 = L1 + math.abs(P1[i] - P0[i]) * dx
   end
   self.L1error = L1
   if self.user_opts.plot then
      util.plot{['code' ]=P     [{{Ng,-Ng},{0,1}}]:table(),
		['exact']=Pexact[{{Ng,-Ng},{0,1}}]:table()}
   end
   if self.user_opts.output then
      local f = io.open(self.user_opts.output, 'w')
      local P0 = P[{{Ng,-Ng},{0,1}}]:table()
      local P1 = P[{{Ng,-Ng},{0,2}}]:table()
      local P2 = P[{{Ng,-Ng},{0,3}}]:table()
      local P3 = P[{{Ng,-Ng},{0,4}}]:table()
      local P4 = P[{{Ng,-Ng},{0,5}}]:table()
      for i=1, self.N do
	 local line = string.format(
	    "%d %+12.8e %+12.8e %+12.8e %+12.8e %+12.8e\n",
	    i, P0[i], P1[i], P2[i], P3[i], P4[i])
	 f:write(line)
      end
      f:close()
   end
end

function NSAccretion:fluid() return 'srmhd' end
function NSAccretion:boundary_conditions() return 'outflow' end

function MaraSimulation:initialize_physics()
   self.Primitive = self.problem:solution()
end

function NSAccretion:solution()
   local Nx = self.simulation.N
   local Ny = self.simulation.N
   local Ng = self.simulation.Ng

   local function pinit(x, y, z)
      local vx = 0.0
      local vy = -y * 0.5
      local Bx = -y / (x^2 + y^2)^0.5
      local By =  x / (x^2 + y^2)^0.5
      return {1, 1, vx, vy, 0, Bx, By, 0}
   end

   local P = array.array{Nx + 2*Ng, Ny + 2*Ng, 8}
   Mara.init_prim(P:buffer(), pinit)
   return P
end

function MaraSimulation:user_work_finish()
   self:checkpoint_write()
end

local user_opts = {
   resolution = 32,
   tmax = 10.0,
   cpi = 0.1,
}
local simulation = MaraSimulation(user_opts)
local problem = NSAccretion(user_opts)

simulation:run(problem)
