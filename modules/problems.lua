
local oo       = require 'class'
local array    = require 'array'
local FishCls  = require 'FishClasses'

local TestProblem     = oo.class('TestProblem')
local TwoStateProblem = oo.class('TwoStateProblem', TestProblem)

local problems = {
   TestProblem = TestProblem,
   TwoStateProblem = TwoStateProblem,

   problems_1d = { },
   problems_2d = { },

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

   SmoothKelvinHelmholtz = oo.class('SmoothKelvinHelmholtz', TestProblem),
   TwoDimensionalImplosion = oo.class('TwoDimensionalImplosion', TestProblem),
   ThrowBlobs = oo.class('ThrowBlobs', TestProblem),
   RelativisticVortex = oo.class('RelativisticVortex', TestProblem),
   JetCavity = oo.class('JetCavity', TestProblem),
   Reconnection = oo.class('Reconnection', TestProblem),
   ShapiroLikeRotator = oo.class('ShapiroLikeRotator', TestProblem),
   MagneticTower = oo.class('MagneticTower', TestProblem)
}
for k,v in pairs(problems) do
   if oo.isclass(v) and oo.issubclass(v, TwoStateProblem) then
      problems.problems_1d[k] = v
   end
   problems.problems_1d['soundwave'] = problems.soundwave
   problems.problems_1d['densitywave'] = problems.densitywave
   problems.problems_1d['collapse1d'] = problems.collapse1d
end


function TestProblem:__init__(user_opts)
   self.user_opts = user_opts or { }
end
function TestProblem:initialize_problem() end
function TestProblem:dynamical_time() return 1.0 end
function TestProblem:finish_time() return 1.0 end
function TestProblem:user_work_iteration() end
function TestProblem:user_work_finish() end
function TestProblem:boundary_conditions() return 'periodic' end
function TestProblem:fluid() return self._fluid or 'nrhyd' end
function TestProblem:adiabatic_index() return 1.4 end


function problems.soundwave:fluid()
   return self.user_opts.self_gravity and 'gravs' or 'nrhyd'
end

function problems.soundwave:initialize_problem()
   self.DomainLength = 1.0
   self.BackgroundDensity = 1.0
   self.FourPiG = 1.0
   self.SoundSpeed  = 0.2
   self.WaveNumber = 4 * math.pi
   self.WaveLength = 2 * math.pi / self.WaveNumber
   self.SoundCrossingTime = self.DomainLength / self.SoundSpeed
   self.JeansLength = 2 * self.SoundSpeed * math.pi / (
      self.FourPiG * self.BackgroundDensity)^0.5
   print("[soundwave] the Jeans length is "..self.JeansLength)
   print("[soundwave] the wave length is "..self.WaveLength)
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
   local d = 0.5
   local P = { }

   P[1] = math.abs(x - 0.5) < d/2 and 1.0 or 1e-6
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
   local Dcen = D[#D/2]

   for i=0,#D-1 do
      if D[i] > Dmax then Dmax = D[i] end
   end
   table.insert(self.max_density, {sim.status.simulation_time, Dcen, Dmax})
end

function problems.collapse1d:user_work_finish()
   local rho0 = 1.0
   local d = 0.5
   local L = 1.0
   local M = d * rho0
   local tau = (M/L)^(-0.5)

   local f = io.open('central-density.dat', 'w')
   for i,val in ipairs(self.max_density) do
      local t = val[1]
      local rho = val[2]
      local real_d = L - (L - d) * math.cosh(t/tau) -- width of the density peak
      local real_rho = M / real_d
      f:write(string.format("%f %f %f\n", t, rho, real_rho))
   end
   f:close()
end

function TwoStateProblem:boundary_conditions()
   return 'outflow'
end

function TwoStateProblem:finish_time()
   return self._finish_time
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

problems.Shocktube1._finish_time = 0.20
problems.Shocktube2._finish_time = 0.05
problems.Shocktube3._finish_time = 0.01
problems.Shocktube4._finish_time = 0.02
problems.Shocktube5._finish_time = 0.02
problems.ContactWave._finish_time = 1.0
problems.SrhdCase1DFIM98._finish_time = 0.1
problems.SrhdCase2DFIM98._finish_time = 0.1
problems.SrhdHardTransverseRAM._finish_time = 0.01

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


function problems.SmoothKelvinHelmholtz:initialize_problem()
   self.vertical_Ek = { }
end

function problems.SmoothKelvinHelmholtz:finish_time()
   return 2.5
end

function problems.SmoothKelvinHelmholtz:solution(x,y,z,t)
   local P0   =  2.5
   local D1   =  1.0
   local D2   =  2.0
   local U1   =  0.5
   local U2   = -0.5
   local dL   =  0.035
   local w0   =  0.010
   local prof = 0.5 * (math.tanh((y - 0.25)/dL) - math.tanh((y - 0.75)/dL))
   local rho  = prof + 1.0
   local vx   = prof - 0.5
   local vy   = w0 * math.sin(4 * math.pi * x)
   return { rho, P0, vx, vy, 0.0 }
end

function problems.SmoothKelvinHelmholtz:user_work_iteration()
   local P = self.simulation.Primitive:vector()
   local t = self.simulation.status.simulation_time
   local E = 0.0
   local n = 0
   for i=0,#P-1,5 do
      local rho = P[i + 0]
      local vy = P[i + 3]
      E = E + 0.5 * rho * vy^2
      n = n + 1
   end

   self.vertical_Ek[t] = E / n
end

function problems.SmoothKelvinHelmholtz:user_work_finish()
   --local util = require 'util'
   --util.pretty_print(self.vertical_Ek)
end


function problems.TwoDimensionalImplosion:finish_time() return 5.0 end
function problems.TwoDimensionalImplosion:boundary_conditions() return 'reflect2d' end
function problems.TwoDimensionalImplosion:solution(x,y,z,t)
   -- http://www.astro.princeton.edu/~jstone/Athena/tests/implode/Implode.html
   if x + y > 0.5 then
      return { 1.000, 1.00, 0.0, 0.0, 0.0 }
   else
      return { 0.125, 0.14, 0.0, 0.0, 0.0 }
   end
end


function problems.ThrowBlobs:initialize_problem()

end

function problems.ThrowBlobs:finish_time()
   return 1.0
end

function problems.ThrowBlobs:boundary_conditions()
   return 'outflow'
end

function problems.ThrowBlobs:solution(x,y,z,t)
   local R0 = 0.75
   local dR = 0.75 * R0
   local b0 = 0.05
   local r1 = ((x - 0.5 - dR/2)^2 + (y - 0.5 - R0 + b0)^2)^0.5
   local r2 = ((x - 0.5 + dR/2)^2 + (y - 0.5 + R0 - b0)^2)^0.5
   local cs = 1.0
   local Ma = 0.5
   local D0 = 1e-1
   local D1 = 1e+1
   local gamma = 1.4
   local rho, vx

   if r1 < R0 then
      rho =  D1
      vx  =  0.5 * Ma * cs
   elseif r2 < R0 then
      rho =  D1
      vx  = -0.5 * Ma * cs
   else
      rho =  D0
      vx  =  0.0
   end

   local pre = D1 * cs^2 / gamma
   return { rho, pre, vx, 0.0, 0.0 }
end

function problems.ThrowBlobs:user_work_iteration()

end

function problems.ThrowBlobs:user_work_finish()

end

function problems.RelativisticVortex:initialize_problem()
end
function problems.RelativisticVortex:finish_time()
   return 2.5
end
function problems.RelativisticVortex:solution(x,y,z,t)
   local Power = math.pow
   local E = math.exp(1)
   local gamma = 1.4
   local r0 = 0.10

   local x = (x - 0.5)
   local y = (y - 0.5)
   local r = (x*x + y*y)^0.5

   local P0   =  0.5 -- pressure at the origin
   local D0   =  1.0 -- density at the origin
   local u0   =  1.0 -- u-phi at the origin

   local P = (D0 - D0*gamma + (D0*(-1 + gamma) + gamma*P0)/
	      Power(E,(gamma*(2*r + r0 - Power(E,(2*r)/r0)*r0)*Power(u0,2))/Power(E,(2*r)/r0)/
		    (4.*(-1 + gamma)*r0)))/gamma

   local u    =  u0 * (r/r0) * math.exp(-r/r0)
   local vf   =  u / (1 + u^2)^0.5
   local vx   = -vf * y/r
   local vy   =  vf * x/r

   return { D0, P, vx, vy, 0.0 }
end
function problems.RelativisticVortex:boundary_conditions() return 'periodic' end
function problems.RelativisticVortex:fluid() return 'srhyd' end

function problems.JetCavity:solution(x,y,z,t)

   local x = (x - 0.5)
   local y = (y - 0.5)
   local r = (x*x + y*y)^0.5

   local vy
   local D
   local P=1.0

   if r > 0.25 then
      D = 100.0
   else
      D = 1.0
   end

   if r < 0.025 then
      vy = 0.999
   end

   return { D, P, 0.0, vy, 0.0, 0, 0, 0 }
end
function problems.JetCavity:boundary_conditions() return 'outflow' end
function problems.JetCavity:fluid() return 'srhyd' end

function problems.Reconnection:initialize_problem(x,y,z,t)
   if self.simulation.cart_rank then
      math.randomseed(self.simulation.cart_rank)
   end
end
function problems.Reconnection:solution(x,y,z,t)
   local x = (x - 0.5)
   local y = (y - 0.5)

   local D=1.0
   local P=1.0
   local Bx

   if math.abs(y) < 0.25 then
      Bx = -1.0
   else
      Bx = 1.0
   end
   local vx = (math.random() - 0.5) * 1e-2
   local vy = (math.random() - 0.5) * 1e-2
   return { D, P, vx, vy, 0.0, Bx, 0.0, 0.0 }
end
function problems.Reconnection:boundary_conditions() return 'periodic' end
function problems.Reconnection:fluid() return 'srmhd' end



function problems.ShapiroLikeRotator:solution(x,y,z,t)
   local x = (x - 0.5)
   local y = (y - 0.5)
   local r = (x*x + y*y)^0.5
   local vx = 0.0
   local vy = 0.0
   local Bx = 0.01
   local By = 0.0
   local r0 = 0.25
   if r < r0 then
      vx = -0.5*y*(1 + (r/r0)^2)
      vy =  0.5*x*(1 + (r/r0)^2)
   end
   return { 1, 1, vx, vy, 0.0, Bx, By, 0.0 }
end
function problems.ShapiroLikeRotator:boundary_conditions() return 'outflow' end
function problems.ShapiroLikeRotator:fluid() return 'srmhd' end



function problems.MagneticTower:solution(x,y,z,t)
   local Power = math.pow
   local E = math.exp(1)
   local gamma = 1.4
   local r0 = 0.15
   local z0 = 0.10
   local x = (x - 0.5)
   local y = (y - 0.5)
   local z = (z - 0.5)
   local r = (x^2 + y^2)^0.5
   local P0   =  0.5 -- pressure at the origin
   local D0   =  1.0 -- density at the origin
   local u0   =  1.0 -- / (1.0 + (z/z0)^2) -- u-phi at the origin
   local P = (D0 - D0*gamma + (D0*(-1 + gamma) + gamma*P0)/
	      Power(E,(gamma*(2*r + r0 - Power(E,(2*r)/r0)*r0)*Power(u0,2))/
		    Power(E,(2*r)/r0)/(4.*(-1 + gamma)*r0)))/gamma
   local u    =  u0 * (r/r0) * math.exp(-r/r0)
   local vf   =  u / (1 + u^2)^0.5
   local vx   = -vf * y/r
   local vy   =  vf * x/r
   local Bx = 0.0
   local By = 0.0
   local Bz = 0.0
   return { D0, P, vx, vy, 0.0, Bx, By, Bz }
end
function problems.MagneticTower:boundary_conditions() return 'outflow' end
function problems.MagneticTower:fluid() return 'srhyd' end



return problems
