
local sim      = require 'simulation'
local oo       = require 'class'
local array    = require 'array'
local fish     = require 'fish'
local fluids   = require 'fluids'
local hdf5     = require 'lua-hdf5.LuaHDF5'
local util     = require 'util'

local FishSimulation = oo.class('FishSimulation', sim.SimulationBase)

function FishSimulation:initialize_solver()
   local opts = self.user_opts
   self.CFL = opts.CFL or 0.8
   self.Ng = 3
   self.N = opts.resolution or 128
   self.dx = 1.0 / self.N

   local FL = self.problem:fluid():upper()
   local descr = fluids.descr_new()
   fluids.descr_setfluid(descr, fluids[FL])
   fluids.descr_setgamma(descr, 1.4)
   fluids.descr_seteos(descr, fluids.EOS_GAMMALAW)

   local RS = ('riemann_'..(self.user_opts.riemann or 'hllc')):upper()
   local RC = (self.user_opts.reconstruction or 'plm'):upper()
   local ST = (self.user_opts.solver or 'godunov'):upper()
   local UP = ({midpoint='midpoint',
		rk3='shuosher_rk3'})[self.user_opts.advance or 'rk3']:upper()
   local BC = self.problem:boundary_conditions():upper()

   local scheme = fish.state_new()
   fish.setparami(scheme, fluids[RS], fish.RIEMANN_SOLVER)
   fish.setparami(scheme, fish[RC], fish.RECONSTRUCTION)
   fish.setparami(scheme, fish[ST], fish.SOLVER_TYPE)
   fish.setparami(scheme, fish[BC], fish.BOUNDARY_CONDITIONS)
   fish.setparami(scheme, fish[UP], fish.TIME_UPDATE)
   fish.setparamd(scheme, 2.0, fish.PLM_THETA)

   fish.grav1d_init(descr, self.N)
   fish.grav1d_setscheme(scheme)

   self.descr = descr
   self.scheme = scheme
end

function FishSimulation:report_configuration()
   local scheme = self.scheme
   local enum = array.vector(1, 'int')
   local cfg = { }

   local FishEnums   = { } -- Register the constants for string lookup later on
   for k,v in pairs(fish) do
      if type(v)=='number' then FishEnums[v]=k end
   end
   local FluidsEnums = { }
   for k,v in pairs(fluids) do
      if type(v)=='number' then FluidsEnums[v]=k end
   end

   for _,k in pairs{'RIEMANN_SOLVER',
		    'RECONSTRUCTION',
		    'SOLVER_TYPE',
		    'BOUNDARY_CONDITIONS',
		    'TIME_UPDATE'} do
      fish.getparami(scheme, enum:buffer(), fish[k])
      local val = FishEnums[enum[0]] or FluidsEnums[enum[0]]
      cfg[k:lower()] = val:lower()
   end

   fluids.descr_getfluid(self.descr, enum:buffer())
   cfg['fluid'] = FluidsEnums[enum[0]]:lower()
   cfg['resolution'] = self.N
   cfg['CFL'] = self.CFL

   print('\t***********************************')
   print('\t*      Solver configuration       *')
   print('\t***********************************')
   util.pretty_print(cfg, '\t+ ')
   print('\t***********************************')
end

function FishSimulation:finalize_solver()
   fish.grav1d_finalize()
   fish.state_del(self.scheme)
   fluids.descr_del(self.descr)
end

function FishSimulation:initialize_physics()
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

function FishSimulation:local_mesh_size()
   return self.N
end

return {FishSimulation=FishSimulation}
