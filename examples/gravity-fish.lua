
local array  = require 'array'
local fish   = require 'fish'
local fluids = require 'fluids'
local hdf5   = require 'lua-hdf5.LuaHDF5'
local util   = require 'util'

local DomainLength = 1.0
local BackgroundDensity = 1.0
local FourPiG = 1.0
local SoundSpeed  = 0.05
local JeansLength = 2 * SoundSpeed * math.pi / (FourPiG * BackgroundDensity)^0.5
local WaveNumber = 8 * math.pi
local WaveLength = 2 * math.pi / WaveNumber
local SoundCrossingTime = DomainLength / SoundSpeed

local N = 512
local Ng = 3
local CFL = 0.5
local dx  = 1.0 / N

util.pretty_print{JeansLength=JeansLength, WaveLength=WaveLength}

local function gravitywave(t)
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
   local w0 = cs * k0 --((cs * k0)^2 - FourPiG * D0)^0.5

   for n=0,#Pvec/5-1 do
      local x  = (n - Ng) * dx

      Pvec[5*n + 0] = D0 + D1 * math.cos(k0*x - w0*t)
      Pvec[5*n + 1] = p0 + p1 * math.cos(k0*x - w0*t)
      Pvec[5*n + 2] = u0 + u1 * math.cos(k0*x - w0*t)
      Pvec[5*n + 3] = 0.0
      Pvec[5*n + 4] = 0.0
   end
   return P, G
end


local P, G = gravitywave(0.0)

fish.grav1d_init(N)
fish.grav1d_mapbuffer(P:buffer(), fluids.PRIMITIVE)
fish.grav1d_mapbuffer(G:buffer(), fluids.GRAVITY)

local t = 0.0
local n = 0

while t < 0.05 * SoundCrossingTime do

   local Amax = fish.grav1d_maxwavespeed()
   local dt = CFL * dx / Amax

   if n % 10 == 0 then
      local fname = string.format('data/gravity-wave-%04d.h5', n)
      local outfile = hdf5.File(fname, 'w')
      outfile['prim'] = P[{{Ng,-Ng},nil}]
      outfile['exact'] = gravitywave(t)[{{Ng,-Ng},nil}]
      outfile:close()
   end

   local start = os.clock()
   fish.grav1d_advance(dt)
   local kzps = 1e-3 * N / (os.clock() - start)

   t = t + dt
   n = n + 1
   if n % 10 == 0 then
      print(string.format("%05d: t=%3.2f dt=%2.1e %4.3fkz/s", n, t, dt, kzps))
   end
end

fish.grav1d_finalize()
