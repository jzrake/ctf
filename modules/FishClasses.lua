
local oo     = require 'class'
local array  = require 'array'
local fish   = require 'fish'
local fluids = require 'fluids'
local util   = require 'util'

local FluidDescriptor = oo.class('FluidDescriptor')
local FluidState      = oo.class('FluidState')
local RiemannSolver   = oo.class('RiemannSolver')
local FishScheme      = oo.class('FishScheme')

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

function FishScheme:__init__(args)
   -- **************************************************************************
   -- args: {table}
   --
   --  + bc: string, boundary conditions [*periodic*|outflow]
   --  + riemann: string, the name of a Riemann solver [hll|*hllc*|exact]
   --  + reconstruction: string [pcm|*plm*|weno5]
   --  + solver: string [*godunov*|spectral]
   --  + advance: string [single|rk2|*rk3*]
   --
   -- **************************************************************************
   local BC = (args.bc or 'periodic'):upper()
   local RS = ('riemann_'..(args.riemann or 'hllc')):upper()
   local RC = (args.reconstruction or 'plm'):upper()
   local ST = (args.solver or 'godunov'):upper()
   local UP = ({ single  = 'single',
		 rk2     = 'tvd_rk2',
		 rk3     = 'shuosher_rk3'
	       })[args.advance or 'rk3']:upper()
   local scheme = fish.state_new()
   fish.setparami(scheme, fluids[RS], fish.RIEMANN_SOLVER)
   fish.setparami(scheme, fish[RC], fish.RECONSTRUCTION)
   fish.setparami(scheme, fish[ST], fish.SOLVER_TYPE)
   fish.setparami(scheme, fish[BC], fish.BOUNDARY_CONDITIONS)
   fish.setparami(scheme, fish[UP], fish.TIME_UPDATE)
   fish.setparamd(scheme, 2.0, fish.PLM_THETA)
   self._c = scheme
end

function FishScheme:__gc__(args)
   fish.state_del(self._c)
end

function FishScheme:report_configuration(extras)
   local scheme = self._c
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

   for k,v in pairs(extras or { }) do
      cfg[k] = v
   end

   cfg['resolution'] = self.N
   cfg['CFL'] = self.CFL

   print('\t***********************************')
   print('\t*      Solver configuration       *')
   print('\t***********************************')
   util.pretty_print(cfg, '\t+ ')
   print('\t***********************************')
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
	   RiemannSolver   = RiemannSolver,
	   FishScheme      = FishScheme}
else
   test1()
   print(debug.getinfo(1).source, ": All tests passed")
end
