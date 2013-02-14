
local sim      = require 'simulation'
local oo       = require 'class'
local array    = require 'array'
local util     = require 'util'
local hdf5     = require 'lua-hdf5.LuaHDF5'
local Mara     = require 'Mara'


local MaraSimulation = oo.class('MaraSimulation', sim.SimulationBase)

function MaraSimulation:initialize_solver()
   local opts = self.user_opts
   self.CFL = opts.CFL or 0.8
   self.Ng = 3
   self.N = opts.resolution or 128
   self.dx = 1.0 / self.N

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
   local N  = self.N
   local Ng = self.Ng
   local dx = self.dx

   local P = array.array{N + 2*Ng, 5}
   local G = array.array{N + 2*Ng, 4}
   local Pvec = P:vector()

   for n=0,#Pvec/5-1 do
      local x  = (n - Ng + 0.5) * dx
      local Pi = self.problem:solution(x,0,0,0)
      Pvec[5*n + 0] = Pi[1]
      Pvec[5*n + 1] = Pi[2]
      Pvec[5*n + 2] = Pi[3]
      Pvec[5*n + 3] = Pi[4]
      Pvec[5*n + 4] = Pi[5]
   end

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

function MaraSimulation:local_mesh_size()
   return self.N
end

return {MaraSimulation=MaraSimulation}
