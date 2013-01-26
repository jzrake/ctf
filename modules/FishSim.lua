
local sim      = require 'simulation'
local oo       = require 'class'
local array    = require 'array'
local fish     = require 'fish'
local fluids   = require 'fluids'
local hdf5     = require 'lua-hdf5.LuaHDF5'
local util     = require 'util'


local FishSimulation = oo.class('FishSimulation', sim.SimulationBase)

function FishSimulation:initialize_behavior()
   local opts = self.user_opts
   local cpi = opts.cpi or 0.1
   local tmax = opts.tmax or 1.0
   local dynamical_time = self.problem:dynamical_time()
   self.behavior.message_cadence = 100
   self.behavior.checkpoint_cadence = cpi * dynamical_time
   self.behavior.max_simulation_time = tmax * dynamical_time
end

function FishSimulation:initialize_solver()
   local opts = self.user_opts
   self.CFL = 0.8
   self.Ng = 3
   self.N = opts.N or 128
   self.dx = 1.0 / self.N

   self.P = array.array{self.N + 2*self.Ng, 5}
   self.G = array.array{self.N + 2*self.Ng, 4}

   local descr = fluids.descr_new()
   fluids.descr_setfluid(descr, fluids.GRAVS)
   fluids.descr_setgamma(descr, 1.4)
   fluids.descr_seteos(descr, fluids.EOS_GAMMALAW)

   local scheme = fish.state_new()
   fish.setparami(scheme, fish.PLM, fish.RECONSTRUCTION)
   fish.setparami(scheme, fluids.RIEMANN_HLLC, fish.RIEMANN_SOLVER)
   fish.setparami(scheme, fish.GODUNOV, fish.SOLVER_TYPE)
   fish.setparamd(scheme, 2.0, fish.PLM_THETA)

   fish.grav1d_init(descr, self.N)
   fish.grav1d_setscheme(scheme)

   self.descr = descr
   self.scheme = scheme
end

function FishSimulation:finalize_solver()
   fish.grav1d_finalize()
   fish.del(self.scheme)
   fluids.descr_del(self.descr)
end

function FishSimulation:initialize_physics()
   local P, G = self.problem:solution(0.0)
   fish.grav1d_mapbuffer(P:buffer(), fluids.PRIMITIVE)
   fish.grav1d_mapbuffer(G:buffer(), fluids.GRAVITY)
   self.Primitive = P
   self.Gravity = G
end

function FishSimulation:set_time_increment()
   local Amax = fish.grav1d_maxwavespeed()
   local dt = self.CFL * self.dx / Amax
   self.status.time_increment = dt
end

function FishSimulation:advance_physics()
   fish.grav1d_advance(self.status.time_increment)
end

function FishSimulation:checkpoint_write()
   local n = self.status.iteration_number
   local t = self.status.simulation_time
   local Ng = self.Ng
   local fname = string.format('data/chkpt.%04d.h5', n)
   local outfile = hdf5.File(fname, 'w')
   outfile['prim' ] = self.Primitive          [{{Ng,-Ng},nil}]
   outfile['grav' ] = self.Gravity            [{{Ng,-Ng},nil}]
   outfile['exact'] = self.problem:solution(t)[{{Ng,-Ng},nil}]
   outfile:close()
end

function FishSimulation:user_work_finish()
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
end

return {FishSimulation=FishSimulation}
