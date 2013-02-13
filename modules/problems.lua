
local oo       = require 'class'
local array    = require 'array'
local FishCls  = require 'FishClasses'

local TestProblem     = oo.class('TestProblem')
local TwoStateProblem = oo.class('TwoStateProblem', TestProblem)

local problems = {
   TestProblem = TestProblem,
   TwoStateProblem = TwoStateProblem,

   soundwave   = oo.class('soundwave'  , TestProblem),
   densitywave = oo.class('densitywave', TestProblem),
   collapse1d  = oo.class('collapse1d' , TestProblem),

   Shocktube1  = oo.class('Shocktube1' , TwoStateProblem),
   Shocktube2  = oo.class('Shocktube2' , TwoStateProblem),
   Shocktube3  = oo.class('Shocktube3' , TwoStateProblem),
   Shocktube4  = oo.class('Shocktube4' , TwoStateProblem),
   Shocktube5  = oo.class('Shocktube5' , TwoStateProblem),
   ContactWave = oo.class('ContactWave', TwoStateProblem),

   SrhdCase1DFIM98       = oo.class('SrhdCase1DFIM98'      , TwoStateProblem),
   SrhdCase2DFIM98       = oo.class('SrhdCase2DFIM98'      , TwoStateProblem),
   SrhdHardTransverseRAM = oo.class('SrhdHardTransverseRAM', TwoStateProblem),
}

function TestProblem:__init__(user_opts)
   self.user_opts = user_opts or { }
   self:initialize_problem()
end
function TestProblem:initialize_problem() end
function TestProblem:dynamical_time() return 1.0 end
function TestProblem:finish_time() return 1.0 end
function TestProblem:user_work_iteration() end
function TestProblem:user_work_finish() end
function TestProblem:boundary_conditions() return 'periodic' end
function TestProblem:fluid() return self._fluid or 'nrhyd' end


function problems.soundwave:fluid()
   return self.user_opts.self_gravity and 'gravs' or 'nrhyd'
end

function problems.soundwave:initialize_problem()
   self.DomainLength = 1.0
   self.BackgroundDensity = 1.0
   self.FourPiG = 1.0
   self.SoundSpeed  = 0.2
   self.WaveNumber = 8 * math.pi
   self.WaveLength = 2 * math.pi / self.WaveNumber
   self.SoundCrossingTime = self.DomainLength / self.SoundSpeed
   self.JeansLength = 2 * self.SoundSpeed * math.pi / (
      self.FourPiG * self.BackgroundDensity)^0.5
   --print("problems.soundwave: the Jeans length is "..self.JeansLength)
end

function problems.soundwave:dynamical_time()
   return self.SoundCrossingTime
end

function problems.soundwave:finish_time()
   return self.SoundCrossingTime
end

function problems.soundwave:solution(x,y,z,t)
   local grav = self.user_opts.self_gravity
   local cs = self.SoundSpeed
   local D0 = self.BackgroundDensity
   local D1 = D0 * 1e-6
   local p0 = D0 * cs^2 / 1.4
   local p1 = D1 * cs^2
   local u0 = 0.0
   local u1 = cs * D1 / D0
   local k0 = self.WaveNumber
   local w0 = grav and ((cs * k0)^2 - self.FourPiG * D0)^0.5 or cs * k0
   local P = { }

   P[1] = D0 + D1 * math.cos(k0*x - w0*t)
   P[2] = p0 + p1 * math.cos(k0*x - w0*t)
   P[3] = u0 + u1 * math.cos(k0*x - w0*t)
   P[4] = 0.0
   P[5] = 0.0

   return P
end

function problems.densitywave:dynamical_time()
   return 1.0
end

function problems.densitywave:finish_time()
   return 1.0
end

function problems.densitywave:solution(x,y,z,t)
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

function problems.collapse1d:fluid()
   return'gravs'
end

function problems.collapse1d:initialize_problem()
   self.max_density = { }
end

function problems.collapse1d:solution(x,y,z,t)
   local p0 = 1e-6

   P[1] = math.abs(x - 0.5) < 0.25 and 1.0 or 1e-6
   P[2] = p0
   P[3] = 0.0
   P[4] = 0.0
   P[5] = 0.0

   return P
end

function problems.collapse1d:user_work_iteration()
   local sim = self.simulation
   local D = sim.Primitive[{nil,{0,1}}]:vector()
   local Dmax = 0.0

   if not self.max_density then self.max_density = { } end

   for i=0,#D-1 do
      if D[i] > Dmax then Dmax = D[i] end
   end
   self.max_density[sim.status.simulation_time] = Dmax
end

function problems.collapse1d:user_work_iteration()
   local f = io.open('stuff.dat', 'w')
   for k,v in pairs(self.max_density) do
      f:write(k..' '..v,'\n')
   end
   f:close()
end

function TwoStateProblem:boundary_conditions()
   return 'outflow'
end

function TwoStateProblem:finish_time()
   return self.problem_finish_time
end

function TwoStateProblem:solution(x,y,z,t)

   if t > 0.0 then
      return self:_exact_solution(x,y,z,t)
   end

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

function TwoStateProblem:_exact_solution(x,y,z,t)
   --
   -- Return the exact solution to the Riemann problem defined by state1 and
   -- state2. Assumes gamma=1.4, and only works when fluid='nrhyd'.
   --
   if self:fluid() ~= 'nrhyd' then
      -- Exact solution to Riemann problem only available for fluid='nrhyd'.
      return {0,0,0,0,0}
   end

   local D = FishCls.FluidDescriptor{gamma=1.4, fluid='nrhyd'}
   local R = FishCls.RiemannSolver()
   local SL = FishCls.FluidState(D)
   local SR = FishCls.FluidState(D)

   for i=1,5 do
      SL.primitive[i-1] = self.state1[i]
      SR.primitive[i-1] = self.state2[i]
   end

   return R:solve(SL, SR, (x-0.5)/t):table()
end

problems.Shocktube1.problem_finish_time = 0.20
problems.Shocktube2.problem_finish_time = 0.05
problems.Shocktube3.problem_finish_time = 0.01
problems.Shocktube4.problem_finish_time = 0.02
problems.Shocktube5.problem_finish_time = 0.02
problems.ContactWave.problem_finish_time = 1.0
problems.SrhdCase1DFIM98.problem_finish_time = 0.1
problems.SrhdCase2DFIM98.problem_finish_time = 0.1
problems.SrhdHardTransverseRAM.problem_finish_time = 0.01

problems.SrhdCase1DFIM98._fluid = 'srhyd'
problems.SrhdCase2DFIM98._fluid = 'srhyd'
problems.SrhdHardTransverseRAM._fluid = 'srhyd'

problems.Shocktube1.state1 = {1.000, 1.000, 0.000, 0.0, 0.0}
problems.Shocktube1.state2 = {0.125, 0.100, 0.000, 0.0, 0.0}
problems.Shocktube2.state1 = {1.000, 0.400,-2.000, 0.0, 0.0}
problems.Shocktube2.state2 = {1.000, 0.400, 2.000, 0.0, 0.0}
problems.Shocktube3.state1 = {1.0, 1e+3, 0.0, 0.0, 0.0}
problems.Shocktube3.state2 = {1.0, 1e-2, 0.0, 0.0, 0.0}
problems.Shocktube4.state1 = {1.0, 1e-2, 0.0, 0.0, 0.0}
problems.Shocktube4.state2 = {1.0, 1e+2, 0.0, 0.0, 0.0}
problems.Shocktube5.state1 = {5.99924, 460.894, 19.59750, 0.0, 0.0}
problems.Shocktube5.state2 = {5.99924,  46.095, -6.19633, 0.0, 0.0}
problems.ContactWave.state1 = {1.0, 1.0, 0.0, 0.7, 0.2}
problems.ContactWave.state2 = {0.1, 1.0, 0.0, 0.7, 0.2}

problems.SrhdCase1DFIM98.state1 = {10.0, 13.30, 0.0, 0.0, 0.0}
problems.SrhdCase1DFIM98.state2 = { 1.0,  1e-6, 0.0, 0.0, 0.0}
problems.SrhdCase2DFIM98.state1 = {1, 1e+3, 0.0, 0.0, 0.0}
problems.SrhdCase2DFIM98.state2 = {1, 1e-2, 0.0, 0.0, 0.0}
problems.SrhdHardTransverseRAM.state1 = {1, 1e+3, 0.0, 0.9, 0.0}
problems.SrhdHardTransverseRAM.state2 = {1, 1e-2, 0.0, 0.9, 0.0}

local p1d = { }
for k,v in pairs(problems) do
   if v ~= TestProblem and v ~= TwoStateProblem then
      p1d[k] = v
   end
end

problems.problems_1d = p1d
return problems
