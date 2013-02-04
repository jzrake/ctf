
local sim      = require 'simulation'
local oo       = require 'class'
local array    = require 'array'
local fish     = require 'fish'
local fluids   = require 'fluids'
local hdf5     = require 'lua-hdf5.LuaHDF5'
local util     = require 'util'
local problems = require 'problems'


local FishEnums   = { } -- Register the constants for string lookup later on
local FluidsEnums = { }
for k,v in pairs(fluids) do if type(v)=='number' then FluidsEnums[v]=k end end
for k,v in pairs(fish)   do if type(v)=='number' then   FishEnums[v]=k end end


local StaticMeshRefinement = oo.class('StaticMeshRefinement', sim.SimulationBase)

function StaticMeshRefinement:initialize_behavior()
   local opts = self.user_opts
   local cpi = opts.cpi or 1.0
   local tmax = opts.tmax or self.problem:finish_time()
   local dynamical_time = self.problem:dynamical_time()
   self.behavior.message_cadence = opts.message_cadence or 10
   self.behavior.checkpoint_cadence = cpi * dynamical_time
   self.behavior.max_simulation_time = tmax * dynamical_time
end

function StaticMeshRefinement:initialize_solver()
   local opts = self.user_opts
   self.CFL = opts.CFL or 0.8
   self.Ng = 3
   self.N = opts.resolution or 128
   self.dx = 1.0 / self.N

   local FL = self.problem:fluid():upper()
   local RS = ('riemann_'..(self.user_opts.riemann or 'hllc')):upper()
   local RC = (self.user_opts.reconstruction or 'plm'):upper()
   local ST = (self.user_opts.solver or 'godunov'):upper()
   local UP = ({single='single',
		midpoint='midpoint',
		rk3='shuosher_rk3'})[self.user_opts.advance or 'rk3']:upper()
   local BC = self.problem:boundary_conditions():upper()

   local descr = fluids.descr_new()
   fluids.descr_setfluid(descr, fluids[FL])
   fluids.descr_setgamma(descr, 1.4)
   fluids.descr_seteos(descr, fluids.EOS_GAMMALAW)

   local blocks = { }
   for i=1,1 do
      local block = fish.block_new()
      fish.block_setdescr(block, descr)
      fish.block_setrank(block, 1)
      fish.block_setsize(block, 0, self.N)
      fish.block_setrange(block, 0, 0.0, 1.0)
      fish.block_setguard(block, self.Ng)
      fish.block_allocate(block)
      blocks[i] = block
   end
   for i=1,1 do
      local BL = blocks[i-1] or blocks[1]
      local BR = blocks[i+1] or blocks[1]
      fish.block_setneighbor(blocks[i], 0, fish.LEFT, BL)
      fish.block_setneighbor(blocks[i], 0, fish.RIGHT, BR)
   end

   local scheme = fish.state_new()
   fish.setparami(scheme, fluids[RS], fish.RIEMANN_SOLVER)
   fish.setparami(scheme, fish[RC], fish.RECONSTRUCTION)
   fish.setparami(scheme, fish[ST], fish.SOLVER_TYPE)
   fish.setparami(scheme, fish[BC], fish.BOUNDARY_CONDITIONS)
   fish.setparami(scheme, fish[UP], fish.TIME_UPDATE)
   fish.setparamd(scheme, 2.0, fish.PLM_THETA)

   self.descr = descr
   self.blocks = blocks
   self.scheme = scheme
   self.primitive = { }
end

function StaticMeshRefinement:report_configuration()
   local scheme = self.scheme
   local enum = array.vector(1, 'long')
   local cfg = { }
   for _,k in pairs{'RIEMANN_SOLVER',
		    'RECONSTRUCTION',
		    'SOLVER_TYPE',
		    'BOUNDARY_CONDITIONS',
		    'TIME_UPDATE'} do
      fish.getparami(scheme, enum:buffer(), fish[k])
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

function StaticMeshRefinement:finalize_solver()
   for _,block in pairs(self.blocks) do
      fish.block_del(self.block)
   end
   fish.state_del(self.scheme)
   fluids.descr_del(self.descr)
end

function StaticMeshRefinement:initialize_physics()
   for _,block in pairs(self.blocks) do
      local P = self.problem:solution(0.0)
      fish.block_mapbuffer(block, P:buffer(), fluids.PRIMITIVE)
      self.primitive[block] = P
   end
end

function StaticMeshRefinement:set_time_increment()
   local Amax = 0.0
   for _,block in pairs(self.blocks) do
      local A = fish.block_maxwavespeed(block)
      if A > Amax then
	 Amax = A
      end
   end
   local dt = self.CFL * self.dx / Amax
   self.status.time_increment = dt
end

function StaticMeshRefinement:advance_physics()

   local dt = self.status.time_increment
   local enum = array.vector(1, 'int')
   local block = self.blocks[1]
   fish.getparami(self.scheme, enum:buffer(), fish.TIME_UPDATE)

   if enum[0] == fish.SINGLE then
      local W0 = array.vector{1.0, 0.0, 1.0}

      for _,block in pairs(self.blocks) do
	 fish.block_fillconserved(block)
      end
      -- ****************************** Step 1 ****************************** --
      for _,block in pairs(self.blocks) do
	 fish.block_timederivative(block, self.scheme)
	 fish.block_evolve(block, W0:buffer(), dt)
      end
      for _,block in pairs(self.blocks) do
	 fish.block_fillguard(block)
      end

   elseif enum[0] == fish.MIDPOINT then
      local W0 = array.vector{1.0, 0.0, 0.5}
      local W1 = array.vector{1.0, 0.0, 1.0}

      for _,block in pairs(self.blocks) do
	 fish.block_fillconserved(block)
      end
      -- ****************************** Step 1 ****************************** --
      for _,block in pairs(self.blocks) do
	 fish.block_timederivative(block, self.scheme)
	 fish.block_evolve(block, W0:buffer(), dt)
      end
      for _,block in pairs(self.blocks) do
	 fish.block_fillguard(block)
      end
      -- ****************************** Step 2 ****************************** --
      for _,block in pairs(self.blocks) do
	 fish.block_timederivative(block, self.scheme)
	 fish.block_evolve(block, W1:buffer(), dt)
      end
      for _,block in pairs(self.blocks) do
	 fish.block_fillguard(block)
      end

   elseif enum[0] == fish.SHUOSHER_RK3 then
      local W0 = array.vector{1.0, 0.0, 1.0}
      local W1 = array.vector{3/4, 1/4, 1/4}
      local W2 = array.vector{1/3, 2/3, 2/3}

      for _,block in pairs(self.blocks) do
	 fish.block_fillconserved(block)
      end
      -- ****************************** Step 1 ****************************** --
      for _,block in pairs(self.blocks) do
	 fish.block_timederivative(block, self.scheme)
	 fish.block_evolve(block, W0:buffer(), dt)
      end
      for _,block in pairs(self.blocks) do
	 fish.block_fillguard(block)
      end
      -- ****************************** Step 2 ****************************** --
      for _,block in pairs(self.blocks) do
	 fish.block_timederivative(block, self.scheme)
	 fish.block_evolve(block, W1:buffer(), dt)
      end
      for _,block in pairs(self.blocks) do
	 fish.block_fillguard(block)
      end
      -- ****************************** Step 3 ****************************** --
      for _,block in pairs(self.blocks) do
	 fish.block_timederivative(block, self.scheme)
	 fish.block_evolve(block, W2:buffer(), dt)
      end
      for _,block in pairs(self.blocks) do
	 fish.block_fillguard(block)
      end
   end
end

function StaticMeshRefinement:checkpoint_write()

end

function StaticMeshRefinement:local_mesh_size()
   return self.N
end

function StaticMeshRefinement:user_work_iteration()
   self.problem:user_work_iteration()
end

function StaticMeshRefinement:user_work_finish()
   local t = self.status.simulation_time
   local Ng = self.Ng
   local P = self.primitive[self.blocks[1]]
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
      local P1 = P[{{Ng,-Ng},{1,2}}]:table()
      local P2 = P[{{Ng,-Ng},{2,3}}]:table()
      local P3 = P[{{Ng,-Ng},{3,4}}]:table()
      local P4 = P[{{Ng,-Ng},{4,5}}]:table()
      for i=1, self.N do
	 local line = string.format(
	    "%d %+12.8e %+12.8e %+12.8e %+12.8e %+12.8e\n",
	    i, P0[i], P1[i], P2[i], P3[i], P4[i])
	 f:write(line)
      end
      f:close()
   end
end

local opts = {plot=true,
	      CFL=0.4,
	      tmax=0.2,
	      solver='godunov',
	      reconstruction='weno5',
	      advance='rk3'}
local sim = StaticMeshRefinement(opts)
local problem = problems.densitywave(opts)
sim:run(problem)
