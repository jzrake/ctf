
local oo    = require 'class'
local array = require 'array'

local problems = { soundwave   = oo.class('soundwave'),
		   densitywave = oo.class('densitywave') }

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

return problems
