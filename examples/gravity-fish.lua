
local array  = require 'array'
local fish   = require 'fish'
local fluids = require 'fluids'
local hdf5   = require 'lua-hdf5.LuaHDF5'

local N = 400
local Ng = 3
local P = array.array{N + 2*Ng, 5}
local G = array.array{N + 2*Ng, 4}
local Pvec = P:vector()

local CFL = 0.5
local dx  = 1.0 / N
local cs  = 1.0

local function gravitywave()
   for n=0,#Pvec/5-1 do
      local x = (n - 50) * dx
      local D0 = 1.0
      local D1 = D0 * 1e-6
      local p0 = D0 * cs^2 / 1.4
      local p1 = D1 * cs^2
      local u1 = cs * D1 / D0
      local k0 = 4 * math.pi

      Pvec[5*n + 0] = D0 + D1 * math.cos(k0 * x)
      Pvec[5*n + 1] = p0 + p1 * math.cos(k0 * x)
      Pvec[5*n + 2] =      u1 * math.cos(k0 * x)
      Pvec[5*n + 3] = 0.0
      Pvec[5*n + 4] = 0.0
   end
end

gravitywave()
fish.grav1d_init(N)
fish.grav1d_mapbuffer(P:buffer(), fluids.PRIMITIVE)

local t = 0.0
local n = 0

while t < 0.3 do

   local Amax = fish.grav1d_maxwavespeed()
   local dt = CFL * dx / Amax

   fish.grav1d_advance(dt)
   t = t + dt
   n = n + 1
   if n % 10 == 0 then
      local fname = string.format('data/euler-%04d.h5', n)
      local outfile = hdf5.File(fname, 'w')
      outfile['prim'] = P[{{Ng,-Ng},nil}]
      outfile:close()
   end
end

fish.grav1d_finalize()
