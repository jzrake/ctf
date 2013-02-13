
local oo     = require 'class'
local array  = require 'array'
local fish   = require 'fish'
local fluids = require 'fluids'

local FluidDescriptor = oo.class('FluidDescriptor')
local FluidState      = oo.class('FluidState')
local RiemannSolver   = oo.class('RiemannSolver')

function FluidDescriptor:__init__(args)
   local args = args or { }
   local descr = fluids.descr_new()
   fluids.descr_setfluid(descr, fluids[(args.fluid or 'nrhyd'):upper()])
   fluids.descr_setgamma(descr, args.gamma or 1.4)
   self._c = descr
end

function FluidDescriptor:__gc__()
   fluids.descr_del(self._c)
end

function FluidDescriptor:get_ncomp(which)
   return fluids.descr_getncomp(self._c, fluids[which:upper()])
end

function FluidDescriptor:fluid()
   local E = { }
   for k,v in pairs(fluids) do
      if type(v) == 'number' then
	 E[v] = k
      end
   end
   local F = array.vector(1, 'int')
   fluids.descr_getfluid(self._c, F:buffer())
   return E[F[0]]:lower()
end

function FluidState:__init__(descr)
   self._d = descr
   self._c = fluids.state_new()
   self._P = array.vector(5)
   fluids.state_setdescr(self._c, descr._c)
   fluids.state_mapbuffer(self._c, self._P:buffer(), fluids.PRIMITIVE)
end

function FluidState:__index__(key)
   return ( -- These are the class read-only attribute
      {descr=self._d,
       primitive=self._P})[key] or oo.getattrib(self, key)  
end

function FluidState:__gc__()
   fluids.state_del(self._c)
end

function RiemannSolver:__init__()
   self._c = fluids.riemn_new()
end

function RiemannSolver:__gc__()
   fluids.riemn_del(self._c)
end

function RiemannSolver:solve(SL, SR, x)
   local S = FluidState(SL.descr)
   fluids.riemn_setsolver (self._c, fluids.RIEMANN_EXACT)
   fluids.riemn_setdim    (self._c, 0)
   fluids.riemn_setstateL (self._c, SL._c)
   fluids.riemn_setstateR (self._c, SR._c)
   fluids.riemn_execute   (self._c)
   fluids.riemn_sample    (self._c, S._c, x)
   return S.primitive
end


local function test1()
   local D = FluidDescriptor()
   local R = RiemannSolver()
   local SL = FluidState(D)
   local SR = FluidState(D)

   SL.primitive[0] = 1.0
   SL.primitive[1] = 1.0
   SR.primitive[0] = 0.125
   SR.primitive[1] = 0.1

   local soln = R:solve(SL, SR, 0.0)
   local real = array.vector{0.42631942817850,
			     0.30313017805065,
			     0.92745262004895, 0, 0}
   for i=0,4 do
      assert(math.abs(soln[i] - real[i]) < 1e-10)
   end
end

if ... then -- if __name__ == "__main__"
   return {FluidDescriptor = FluidDescriptor,
	   FluidState      = FluidState,
	   RiemannSolver   = RiemannSolver}
else
   test1()
   print(debug.getinfo(1).source, ": All tests passed")
end
