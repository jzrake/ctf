
local sim      = require 'simulation'
local oo       = require 'class'
local array    = require 'array'
local fish     = require 'fish'
local fluids   = require 'fluids'
local hdf5     = require 'lua-hdf5.LuaHDF5'
local util     = require 'util'
local mesh     = require 'mesh'
local problems = require 'problems'


local densitywave = oo.class('densitywave', problems.TestProblem)
local ST1 = oo.class('ST1', problems.Shocktube1)
function densitywave:dynamical_time()
   return 1.0
end
function densitywave:finish_time()
   return 1.0
end
function densitywave:solution(x, y, z, t)
   local sim = self.simulation
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

function ST1:solution(x, y, z, t)
   local sim = self.simulation
   local L = self.state1
   local R = self.state2
   local P = { }
   P[1] = x < 0.5 and L[1] or R[1]
   P[2] = x < 0.5 and L[2] or R[2]
   P[3] = x < 0.5 and L[3] or R[3]
   P[4] = x < 0.5 and L[4] or R[4]
   P[5] = x < 0.5 and L[5] or R[5]
   return P
end

local FishEnums   = { } -- Register the constants for string lookup later on
for k,v in pairs(fish) do
   if type(v)=='number' then FishEnums[v]=k end
end
local FluidsEnums = { }
for k,v in pairs(fluids) do
   if type(v)=='number' then FluidsEnums[v]=k end
end

local StaticMeshRefinement = oo.class('StaticMeshRefinement', sim.SimulationBase)
      
function StaticMeshRefinement:initialize_behavior()
   local opts = self.user_opts
   local cpi = opts.cpi or 1.0
   local tmax = opts.tmax or self.problem:finish_time()
   local dynamical_time = self.problem:dynamical_time()
   self.behavior.message_cadence = opts.message_cadence or 1
   self.behavior.checkpoint_cadence = cpi * dynamical_time
   self.behavior.max_simulation_time = tmax * dynamical_time
end


local function TiledUniformLevelMesh(args)
   -- **************************************************************************
   --
   -- Create a mesh for which only the finest level is active, built from blocks
   -- of size `N` to depth `level`
   --
   -- args: {table}
   --
   --  + N: number, size of block
   --  + Ng: number, guard zones
   --  + level: number, depth of the grid
   --  + periodic: boolean, true by default
   --
   -- **************************************************************************
   --
   local mesh = mesh.Block { size={args.N},
			     guard=args.guard,
			     dummy=args.level ~= 0 }

   local function expand_block(b, level)
      if level == 0 then return end
      b:add_child_block(0, {dummy=level ~= 1})
      b:add_child_block(1, {dummy=level ~= 1})
      expand_block(b[0], level - 1)
      expand_block(b[1], level - 1)
   end
   expand_block(mesh, args.level)

   local left_most, right_most

   for b in mesh:walk() do
      if not b:neighbor_block(0, 'L') then
	 left_most = b
      end
      if not b:neighbor_block(0, 'R') then
	 right_most = b
      end
   end

   if not args.periodic then
      left_most :set_boundary_block(0, 'L', left_most)
      right_most:set_boundary_block(0, 'R', right_most)
   else
      left_most :set_boundary_block(0, 'L', right_most)
      right_most:set_boundary_block(0, 'R', left_most)
   end

   return mesh
end


function StaticMeshRefinement:initialize_solver()
   local opts = self.user_opts
   self.CFL = opts.CFL or 0.8
   self.Ng = 3
   self.N = opts.resolution or 16

   local BC = self.problem:boundary_conditions():upper()
   local FL = self.problem:fluid():upper()
   local RS = ('riemann_'..(self.user_opts.riemann or 'hllc')):upper()
   local RC = (self.user_opts.reconstruction or 'plm'):upper()
   local ST = (self.user_opts.solver or 'godunov'):upper()
   local UP = ({ single  = 'single',
		 rk2     = 'tvd_rk2',
		 rk3     = 'shuosher_rk3'
	       })[self.user_opts.advance or 'rk3']:upper()

   local mesh = TiledUniformLevelMesh{ N=self.N,
				       level=3,
				       guard=self.Ng,
				       periodic=true }

   mesh[0][1][0]:add_child_block(0):add_child_block(0)
   mesh[0][1][1]:add_child_block(1):add_child_block(1)

   local scheme = fish.state_new()
   fish.setparami(scheme, fluids[RS], fish.RIEMANN_SOLVER)
   fish.setparami(scheme, fish[RC], fish.RECONSTRUCTION)
   fish.setparami(scheme, fish[ST], fish.SOLVER_TYPE)
   fish.setparami(scheme, fish[BC], fish.BOUNDARY_CONDITIONS)
   fish.setparami(scheme, fish[UP], fish.TIME_UPDATE)
   fish.setparamd(scheme, 2.0, fish.PLM_THETA)

   self.mesh   = mesh
   self.scheme = scheme
end

function StaticMeshRefinement:report_configuration()
   local scheme = self.scheme
   local enum = array.vector(1, 'int')
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

   cfg['fluid'] = self.mesh.descr:fluid()
   cfg['resolution'] = self.N
   cfg['CFL'] = self.CFL

   print('\t***********************************')
   print('\t*      Solver configuration       *')
   print('\t***********************************')
   util.pretty_print(cfg, '\t+ ')
   print('\t***********************************')
end

function StaticMeshRefinement:finalize_solver()
   fish.state_del(self.scheme)
end

function StaticMeshRefinement:initialize_physics()
   local function primitive_init(x,y,z)
      return self.problem:solution(x,y,z,0)
   end
   for b in self.mesh:walk() do
      b:map(primitive_init)
   end
end

function StaticMeshRefinement:set_time_increment()
   local Amax = self.mesh:max_wavespeed()
   local Dmin = self.mesh:grid_spacing{ mode='smallest' }
   local dt = self.CFL * Dmin / Amax
   self.status.time_increment = dt
end

function StaticMeshRefinement:advance_physics()
   local dt = self.status.time_increment
   local enum = array.vector(1, 'int')
   fish.getparami(self.scheme, enum:buffer(), fish.TIME_UPDATE)

   if enum[0] == fish.SINGLE then
      local W0 = array.vector{ 0.0, 1.0 }

      for block in self.mesh:walk() do
	 block:fill_conserved()
      end
      -- ****************************** Step 1 ****************************** --
      for block in self.mesh:walk() do
	 block:time_derivative(self.scheme)
	 block:evolve(W0, dt)
      end
      self.mesh:fill()
      for block in self.mesh:walk() do
	 block:fill_guard(block)
      end

   elseif enum[0] == fish.TVD_RK2 then
      local W0 = array.vector{ 0.0, 1.0 }
      local W1 = array.vector{ 0.5, 0.5 }

      for block in self.mesh:walk() do
	 block:fill_conserved()
      end
      -- ****************************** Step 1 ****************************** --
      for block in self.mesh:walk() do
	 block:time_derivative(self.scheme)
	 block:evolve(W0, dt)
      end
      self.mesh:fill()
      for block in self.mesh:walk() do
	 block:fill_guard(block)
      end
      -- ****************************** Step 2 ****************************** --
      for block in self.mesh:walk() do
	 block:time_derivative(self.scheme)
	 block:evolve(W1, dt)
      end
      self.mesh:fill()
      for block in self.mesh:walk() do
	 block:fill_guard(block)
      end

   elseif enum[0] == fish.SHUOSHER_RK3 then
      local W0 = array.vector{ 0.0, 1.0 }
      local W1 = array.vector{ 3/4, 1/4 }
      local W2 = array.vector{ 1/3, 2/3 }

      for block in self.mesh:walk() do
	 block:fill_conserved(block)
      end
      -- ****************************** Step 1 ****************************** --
      for block in self.mesh:walk() do
	 block:time_derivative(self.scheme)
	 block:evolve(W0, dt)
      end
      self.mesh:fill()
      for block in self.mesh:walk() do
	 block:fill_guard(block)
      end
      -- ****************************** Step 2 ****************************** --
      for block in self.mesh:walk() do
	 block:time_derivative(self.scheme)
	 block:evolve(W1, dt)
      end
      self.mesh:fill()
      for block in self.mesh:walk() do
	 block:fill_guard(block)
      end
      -- ****************************** Step 3 ****************************** --
      for block in self.mesh:walk() do
	 block:time_derivative(self.scheme)
	 block:evolve(W2, dt)
      end
      self.mesh:fill()
      for block in self.mesh:walk() do
	 block:fill_guard(block)
      end
   end
end

function StaticMeshRefinement:checkpoint_write()

end

function StaticMeshRefinement:local_mesh_size()
   return self.mesh:total_states{recurse=true, mode='interior'}
end

function StaticMeshRefinement:user_work_iteration()
   self.problem:user_work_iteration()
end

function StaticMeshRefinement:user_work_finish()
   local t = self.status.simulation_time

   local levels = { }
   local blocks = { }
   local all = { }
   for b in self.mesh:walk() do
      local D = b:level()
      levels[D] = b:table(0, levels[D])
      blocks[b] = b:table(0)
      b:table(0, all)
   end

   if self.user_opts.plot then
      --util.plot({all=all}, {ls='w p', output=nil})
      util.plot(levels, {ls='w p', output=nil})
      --util.plot(blocks, {ls='w p', output=nil})
   end
end

local opts = {plot=true,
	      resolution=16,
	      CFL=0.8,
	      tmax=0.1,
	      solver='spectral',
	      reconstruction='weno5',
	      advance='rk3'}

local sim = StaticMeshRefinement(opts)
--local problem = ST1(opts)
local problem = densitywave(opts)
sim:run(problem)
