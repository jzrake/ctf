
local oo    = require 'class'
local array = require 'array'

local TestProblem = oo.class('TestProblem')
local TwoStateProblem = oo.class('TwoStateProblem', TestProblem)

local problems = {
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
   self.user_opts = user_opts
end
function TestProblem:dynamical_time() return 1.0 end
function TestProblem:user_work_iteration() end
function TestProblem:user_work_finish() end
function TestProblem:boundary_conditions()
   return 'periodic'
end

local DomainLength = 1.0
local BackgroundDensity = 1.0
local FourPiG = 1.0
local SoundSpeed  = 0.2
local JeansLength = 2 * SoundSpeed * math.pi / (FourPiG * BackgroundDensity)^0.5
local WaveNumber = 8 * math.pi
local WaveLength = 2 * math.pi / WaveNumber
local SoundCrossingTime = DomainLength / SoundSpeed

function problems.soundwave:dynamical_time()
   return SoundCrossingTime
end
function problems.soundwave:solution(t)
   local sim = self.simulation
   local N = sim.N
   local Ng = sim.Ng
   local dx = sim.dx

   local P = array.array{N + 2*Ng, 5}
   local G = array.array{N + 2*Ng, 4}
   local Pvec = P:vector()

   local cs = SoundSpeed
   local D0 = BackgroundDensity
   local D1 = D0 * 1e-6
   local p0 = D0 * cs^2 / 1.4
   local p1 = D1 * cs^2
   local u0 = 0.0
   local u1 = cs * D1 / D0
   local k0 = WaveNumber
   local w0 = ((cs * k0)^2 - FourPiG * D0)^0.5

   for n=0,#Pvec/5-1 do
      local x = (n - sim.Ng) * dx
      Pvec[5*n + 0] = D0 + D1 * math.cos(k0*x - w0*t)
      Pvec[5*n + 1] = p0 + p1 * math.cos(k0*x - w0*t)
      Pvec[5*n + 2] = u0 + u1 * math.cos(k0*x - w0*t)
      Pvec[5*n + 3] = 0.0
      Pvec[5*n + 4] = 0.0
   end
   return P, G
end

function problems.densitywave:dynamical_time()
   return SoundCrossingTime
end
function problems.densitywave:solution(t)
   local sim = self.simulation
   local N = sim.N
   local Ng = sim.Ng
   local dx = sim.dx

   local P = array.array{N + 2*Ng, 5}
   local G = array.array{N + 2*Ng, 4}
   local Pvec = P:vector()

   local cs = SoundSpeed
   local u0 = cs -- Mach 1 density wave
   local D0 = BackgroundDensity
   local D1 = D0 * 1e-1
   local p0 = D0 * cs^2 / 1.4
   local k0 = WaveNumber

   for n=0,#Pvec/5-1 do
      local x  = (n - Ng) * dx
      Pvec[5*n + 0] = D0 + D1 * math.cos(k0 * (x - u0 * t))
      Pvec[5*n + 1] = p0
      Pvec[5*n + 2] = u0
      Pvec[5*n + 3] = 0.0
      Pvec[5*n + 4] = 0.0
   end
   return P, G
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

function problems.collapse1d:solution(t)
   local sim = self.simulation
   local N = sim.N
   local Ng = sim.Ng
   local dx = sim.dx
   local p0 = 1e-6

   local P = array.array{N + 2*Ng, 5}
   local G = array.array{N + 2*Ng, 4}
   local Pvec = P:vector()

   for n=0,#Pvec/5-1 do
      local x  = (n - Ng) * dx
      Pvec[5*n + 0] = math.abs(x - 0.5) < 0.25 and 1.0 or 1e-6
      Pvec[5*n + 1] = p0
      Pvec[5*n + 2] = 0.0
      Pvec[5*n + 3] = 0.0
      Pvec[5*n + 4] = 0.0
   end
   return P, G
end



function TwoStateProblem:boundary_conditions() 
   return 'outflow'
end
function TwoStateProblem:solution(t)
   local sim = self.simulation
   local N = sim.N
   local Ng = sim.Ng
   local dx = sim.dx
   local p0 = 1e-6

   local L = self.state1
   local R = self.state2

   local P = array.array{N + 2*Ng, 5}
   local G = array.array{N + 2*Ng, 4}
   local Pvec = P:vector()

   for n=0,#Pvec/5-1 do
      local x  = (n - Ng) * dx
      Pvec[5*n + 0] = x < 0.5 and L[1] or R[1]
      Pvec[5*n + 1] = x < 0.5 and L[2] or R[2]
      Pvec[5*n + 2] = x < 0.5 and L[3] or R[3]
      Pvec[5*n + 3] = x < 0.5 and L[3] or R[4]
      Pvec[5*n + 4] = x < 0.5 and L[4] or R[5]
   end
   return P, G
end

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

return problems
