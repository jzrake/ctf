
local sim      = require 'simulation'
local oo       = require 'class'
local array    = require 'array'
local fish     = require 'fish'
local fluids   = require 'fluids'
local util     = require 'util'
local mesh     = require 'mesh'
local problems = require 'problems'
local FishCls  = require 'FishClasses'


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
   --  + bc: string (boundary conditions), ['periodic', 'outflow']
   --
   -- **************************************************************************
   --
   local descr = FishCls.FluidDescriptor{fluid='nrhyd'}
   local root = mesh.Block { descr=descr,
			     size={args.N},
			     guard=args.guard,
			     dummy=args.level ~= 0 }

   local function expand_block(b, level)
      if level == 0 then return end
      b:add_child_block(0, {dummy=level ~= 1})
      b:add_child_block(1, {dummy=level ~= 1})
      expand_block(b[0], level - 1)
      expand_block(b[1], level - 1)
   end
   expand_block(root, args.level)

   local left_most, right_most

   for b in root:walk() do
      if not b:neighbor_block(0, 'L') then
	 left_most = b
      end
      if not b:neighbor_block(0, 'R') then
	 right_most = b
      end
   end

   if args.bc == 'outflow' then
      left_most :set_boundary_flag(0, 'L', 'outflow')
      right_most:set_boundary_flag(0, 'R', 'outflow')
   elseif args.bc == 'periodic' then
      left_most :set_boundary_block(0, 'L', right_most)
      right_most:set_boundary_block(0, 'R', left_most)
   else
      error("must give bc='periodic or bc='outflow', got "..tostring(bc))
   end

   return root
end


local StaticMeshRefinement = oo.class('StaticMeshRefinement', sim.SimulationBase)

function StaticMeshRefinement:initialize_behavior()
   local opts = self.user_opts
   local cpi = opts.cpi or 1.0
   local tmax = opts.tmax or self.problem:finish_time()
   self.behavior.message_cadence = opts.message_cadence or 1
   self.behavior.checkpoint_cadence = cpi
   self.behavior.max_simulation_time = tmax
end

function StaticMeshRefinement:initialize_solver()
   local opts = self.user_opts
   self.CFL = opts.CFL or 0.8
   self.Ng = 3
   self.N = opts.resolution or 16
   local scheme = FishCls.FishScheme { bc=self.problem:boundary_conditions(),
				       riemann=self.user_opts.riemann,
				       reconstruction=self.user_opts.reconstruction,
				       solver=self.user_opts.solver,
				       advance=self.user_opts.advance }
   local mesh = TiledUniformLevelMesh { N=self.N,
					level=3,
					guard=self.Ng,
					bc=self.problem:boundary_conditions() }
   self.mesh   = mesh
   self.scheme = scheme
end

function StaticMeshRefinement:report_configuration()
   self.scheme:report_configuration()
end

function StaticMeshRefinement:finalize_solver()

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
   fish.getparami(self.scheme._c, enum:buffer(), fish.TIME_UPDATE)

   local function step(w0, w1)
      local W = array.vector{ w0, w1 }
      for block in self.mesh:walk() do
	 block:time_derivative(self.scheme._c)
	 block:source_terms()
	 if self.mesh.descr:fluid() == 'gravs' then
	    block:solve_poisson()
	 end
	 block:evolve(W, dt)
      end
      self.mesh:fill()
      for block in self.mesh:walk() do
	 block:fill_guard(block)
      end
   end

   for block in self.mesh:walk() do
      block:fill_conserved()
   end

   if enum[0] == fish.SINGLE then
      step(0.0, 1.0)

   elseif enum[0] == fish.TVD_RK2 then
      step(0.0, 1.0)
      step(0.5, 0.5)

   elseif enum[0] == fish.SHUOSHER_RK3 then
      step(0.0, 1.0)
      step(3/4, 1/4)
      step(1/3, 2/3)
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
local problem = problems.densitywave(opts)
--local problem = problems.Shocktube1(opts)
sim:run(problem)
