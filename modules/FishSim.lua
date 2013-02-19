
local sim      = require 'simulation'
local oo       = require 'class'
local array    = require 'array'
local fish     = require 'fish'
local fluids   = require 'fluids'
local util     = require 'util'
local mesh     = require 'mesh'
local FishCls  = require 'FishClasses'

local FishSimulation = oo.class('FishSimulation', sim.SimulationBase)


function FishSimulation:initialize_solver()
   local opts = self.user_opts
   self.CFL = opts.CFL or 0.8
   self.Ng = 3
   self.N = opts.resolution or 128
   local bcflag = self.problem:boundary_conditions()
   local descr = FishCls.FluidDescriptor{ fluid=self.problem:fluid() }
   local scheme = FishCls.FishScheme { bc=bcflag,
				       riemann=self.user_opts.riemann,
				       reconstruction=self.user_opts.reconstruction,
				       solver=self.user_opts.solver,
				       advance=self.user_opts.advance }
   local grid = mesh.Block { descr=descr,
			     size={self.N},
			     guard=self.Ng }

   if bcflag == 'periodic' then
      grid:set_boundary_block(0, 'L', grid)
      grid:set_boundary_block(0, 'R', grid)
   else
      grid:set_boundary_flag(0, 'L', bcflag)
      grid:set_boundary_flag(0, 'R', bcflag)
   end

   self.grid   = grid
   self.scheme = scheme
   self.descr  = descr
end

function FishSimulation:report_configuration()
   self.scheme:report_configuration{ fluid=self.descr:fluid() }
end

function FishSimulation:initialize_physics()
   local function primitive_init(x,y,z)
      return self.problem:solution(x,y,z,0)
   end
   for b in self.grid:walk() do
      b:map(primitive_init)
   end
end

function FishSimulation:set_time_increment()
   local Amax = self.grid:max_wavespeed()
   local Dmin = self.grid:grid_spacing{ mode='smallest' }
   local dt = self.CFL * Dmin / Amax
   self.status.time_increment = dt
end

function FishSimulation:advance_physics()
   local dt = self.status.time_increment
   local enum = array.vector(1, 'int')
   local fid = self.descr:fluid()
   fish.getparami(self.scheme._c, enum:buffer(), fish.TIME_UPDATE)

   local function step(w0, w1)
      local W = array.vector{ w0, w1 }
      for block in self.grid:walk() do
	 if fid == 'gravs' or
	    fid == 'gravp' then
	    block:solve_poisson()
	 end
	 block:time_derivative(self.scheme._c)
	 block:source_terms()
	 block:evolve(W, dt)
      end
      self.grid:fill()
      for block in self.grid:walk() do
	 block:fill_guard(block)
      end
   end

   for block in self.grid:walk() do
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

function FishSimulation:local_mesh_size()
   return self.grid:total_states{recurse=true, mode='interior'}
end

return {FishSimulation=FishSimulation}
