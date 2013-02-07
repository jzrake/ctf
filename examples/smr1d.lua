
local sim      = require 'simulation'
local oo       = require 'class'
local array    = require 'array'
local fish     = require 'fish'
local fluids   = require 'fluids'
local hdf5     = require 'lua-hdf5.LuaHDF5'
local util     = require 'util'
local problems = require 'problems'


local densitywave = oo.class('densitywave', problems.TestProblem)
function densitywave:dynamical_time()
   return 1.0
end
function densitywave:finish_time()
   return 1.0
end
function densitywave:solution(x, y, z, t)
   local sim = self.simulation
   local N = sim.N
   local Ng = sim.Ng

   local cs = 1.0
   local u0 = cs -- Mach 1 density wave
   local D0 = 1.0
   local D1 = D0 * 1e-1
   local p0 = D0 * cs^2 / 1.4
   local k0 = 2 * math.pi

   local P = { }
   P[1] = D0 + D1 * math.cos(k0 * (x - u0 * t))
   P[2] = p0
   P[3] = u0
   P[4] = 0.0
   P[5] = 0.0
   return P
end



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
   self.N = opts.resolution or 16

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

   local num_blocks = 8
   local X0 = 0.0
   local X1 = 1.0
   local dX = (X1 - X0) / num_blocks

   local blocks = { }
   for i=0,num_blocks-1 do
      local block = fish.block_new()
      fish.block_setdescr(block, descr)
      fish.block_setrank(block, 1)
      fish.block_setsize(block, 0, self.N)
      fish.block_setrange(block, 0, X0 + (i+0)*dX, X0 + (i+1)*dX)
      fish.block_setguard(block, self.Ng)
      fish.block_allocate(block)
      table.insert(blocks, block)
   end
   for i,block in ipairs(blocks) do
      local BL = blocks[i-1] or blocks[num_blocks]
      local BR = blocks[i+1] or blocks[1]
      fish.block_setneighbor(blocks[i], 0, fish.LEFT, BL)
      fish.block_setneighbor(blocks[i], 0, fish.RIGHT, BR)
   end

   
   local blockL = fish.block_new()
   local blockR = fish.block_new()


   for _,block in pairs{blockL, blockR} do
      fish.block_setdescr(block, descr)
      fish.block_setrank(block, 1)
      fish.block_setsize(block, 0, self.N)
      fish.block_setguard(block, self.Ng)
      fish.block_allocate(block)
   end

   fish.block_setneighbor(blockL, 0, fish.RIGHT, blockR)
   fish.block_setneighbor(blockR, 0, fish.LEFT, blockL)
   fish.block_setchild(blocks[3], 0, blockL)
   fish.block_setchild(blocks[3], 1, blockR)

   table.insert(blocks, blockL)
   table.insert(blocks, blockR)

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
      local Nx = fish.block_getsize(block, 0)
      local Ng = fish.block_getguard(block)
      local P = array.array{Nx + 2*Ng, 5}
      local Pvec = P:vector()

      for i=0,Nx+2*Ng-1 do
	 local x = fish.block_positionatindex(block, 0, i)
	 local Pi = self.problem:solution(x, 0.0, 0.0, 0.0)
	 Pvec[5*i + 0] = Pi[1]
	 Pvec[5*i + 1] = Pi[2]
	 Pvec[5*i + 2] = Pi[3]
	 Pvec[5*i + 3] = Pi[4]
	 Pvec[5*i + 4] = Pi[5]
      end

      fish.block_mapbuffer(block, P:buffer(), fluids.PRIMITIVE)
      self.primitive[block] = P
   end
end

function StaticMeshRefinement:set_time_increment()
   local Amax = 0.0
   local Dmin = math.huge

   for _,block in pairs(self.blocks) do
      local D = fish.block_gridspacing(block, 0)
      local A = fish.block_maxwavespeed(block)

      if A > Amax then Amax = A end
      if D < Dmin then Dmin = D end
   end

   local dt = self.CFL * Dmin / Amax
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
   return #self.blocks * self.N
end

function StaticMeshRefinement:user_work_iteration()
   self.problem:user_work_iteration()
end

function StaticMeshRefinement:user_work_finish()
   local t = self.status.simulation_time

   local code_data = { }
   local exac_data = { }

   for i,block in pairs(self.blocks) do
      local Nx = fish.block_getsize(block, 0)
      local Ng = fish.block_getguard(block)

      if i ~= 3 then

	 for i=Ng,Nx+Ng-1 do
	    local x = fish.block_positionatindex(block, 0, i)
	    
	    local P0 = self.primitive[block]:vector()
	    local P1 = self.problem:solution(x, 0.0, 0.0, t)
	    
	    code_data[x] = P0[5*i + 0]
	    exac_data[x] = P1[1]
	 end

      end
   end

   if self.user_opts.plot then
      util.plot({['code' ]=code_data}, {ls='w p', output=nil})
		--['exact']=exac_data}
   end
end

local opts = {plot=true,
	      CFL=0.8,
	      tmax=0.1,
	      solver='godunov',
	      reconstruction='plm',
	      advance='rk3'}
local sim = StaticMeshRefinement(opts)
local problem = densitywave(opts)
sim:run(problem)
