
local sim      = require 'simulation'
local oo       = require 'class'
local array    = require 'array'
local fish     = require 'fish'
local fluids   = require 'fluids'
local hdf5     = require 'lua-hdf5.LuaHDF5'
local util     = require 'util'


local FishEnums   = { } -- Register the constants for string lookup later on
local FluidsEnums = { }
for k,v in pairs(fluids) do if type(v)=='number' then FluidsEnums[v]=k end end
for k,v in pairs(fish)   do if type(v)=='number' then   FishEnums[v]=k end end


local FishSimulation = oo.class('FishSimulation', sim.SimulationBase)

function FishSimulation:initialize_behavior()
   local opts = self.user_opts
   local cpi = opts.cpi or 1.0
   local tmax = opts.tmax or 1.0
   local dynamical_time = self.problem:dynamical_time()
   self.behavior.message_cadence = 10
   self.behavior.checkpoint_cadence = cpi * dynamical_time
   self.behavior.max_simulation_time = tmax * dynamical_time
end

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
   local enum = array.vector(1, 'long')
   local cfg = { }
   for _,k in pairs{'RIEMANN_SOLVER',
		    'RECONSTRUCTION',
		    'SOLVER_TYPE',
		    'BOUNDARY_CONDITIONS',
		    'TIME_UPDATE'} do
      fish.getparami(scheme, enum:pointer(), fish[k])
      local val = FishEnums[enum[0]] or FluidsEnums[enum[0]]
      cfg[k:lower()] = val:lower()
   end

   local enum = array.vector(1, 'int')
   fluids.descr_getfluid(self.descr, enum:pointer())

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

function FishSimulation:user_work_iteration()
   self.problem:user_work_iteration()
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
