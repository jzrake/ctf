
local sim      = require 'simulation'
local oo       = require 'class'
local array    = require 'array'
local util     = require 'util'
local hdf5     = require 'lua-hdf5.LuaHDF5'
local Mara     = require 'Mara'


local MaraSimulation = oo.class('MaraSimulation', sim.SimulationBase)

function MaraSimulation:initialize_behavior()
   local opts = self.user_opts
   local cpi = opts.cpi or 1.0
   local tmax = opts.tmax or 1.0
   local dynamical_time = self.problem:dynamical_time()
   self.behavior.message_cadence = 10
   self.behavior.checkpoint_cadence = cpi * dynamical_time
   self.behavior.max_simulation_time = tmax * dynamical_time
end

function MaraSimulation:initialize_solver()
   local opts = self.user_opts
   self.CFL = opts.CFL or 0.8
   self.Ng = 3
   self.N = opts.N or 128
   self.dx = 1.0 / self.N

   Mara.start()
   Mara.set_fluid('euler')
   Mara.set_advance('rk3')
   Mara.set_godunov('weno-split')
   Mara.set_boundary(self.problem:boundary_conditions())
   Mara.set_riemann('hllc')
   Mara.config_solver({theta=2.0}, true)

   local prim_names = Mara.fluid.GetPrimNames()
   local Nq = #prim_names

   Mara.set_domain({0}, {1}, {self.N}, Nq, self.Ng)
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
   Mara.advance(P, dt)
end

function MaraSimulation:checkpoint_write()
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
end

return {MaraSimulation=MaraSimulation}
