
local array  = require 'array'
local fish   = require 'fish'
local fluids = require 'fluids'
local hdf5   = require 'lua-hdf5.LuaHDF5'
local util   = require 'util'

local OutputCadence = 10
local MessageCadence = 10
local DomainLength = 1.0
local BackgroundDensity = 1.0
local FourPiG = 1.0
local SoundSpeed  = 0.2
local JeansLength = 2 * SoundSpeed * math.pi / (FourPiG * BackgroundDensity)^0.5
local WaveNumber = 8 * math.pi
local WaveLength = 2 * math.pi / WaveNumber
local SoundCrossingTime = DomainLength / SoundSpeed
local ErrorTable = { }

local Ng = 3
local CFL = 0.8
local N, dx

local function soundwave(t)
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
      local x = (n - Ng) * dx
      Pvec[5*n + 0] = D0 + D1 * math.cos(k0*x - w0*t)
      Pvec[5*n + 1] = p0 + p1 * math.cos(k0*x - w0*t)
      Pvec[5*n + 2] = u0 + u1 * math.cos(k0*x - w0*t)
      Pvec[5*n + 3] = 0.0
      Pvec[5*n + 4] = 0.0
   end
   return P, G
end

local function densitywave(t)
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

local solution = soundwave
util.pretty_print{JeansLength=JeansLength, WaveLength=WaveLength}

local function run(NumZones)

   N = NumZones
   dx = 1.0 / N

   local t = 0.0
   local n = 0
   local P, G = solution(0.0)

   local descr = fluids.descr_new()
   fluids.descr_setfluid(descr, fluids.GRAVS)
   fluids.descr_setgamma(descr, 1.4)
   fluids.descr_seteos(descr, fluids.EOS_GAMMALAW)

   local scheme = fish.state_new()
   fish.setparami(scheme, fish.PLM, fish.RECONSTRUCTION)
   fish.setparami(scheme, fluids.RIEMANN_HLLC, fish.RIEMANN_SOLVER)
   fish.setparami(scheme, fish.GODUNOV, fish.SOLVER_TYPE)
   fish.setparamd(scheme, 2.0, fish.PLM_THETA)

   fish.grav1d_init(descr, N)
   fish.grav1d_setscheme(scheme)
   fish.grav1d_mapbuffer(P:buffer(), fluids.PRIMITIVE)
   fish.grav1d_mapbuffer(G:buffer(), fluids.GRAVITY)

   while t < SoundCrossingTime do

      local Amax = fish.grav1d_maxwavespeed()
      local dt = CFL * dx / Amax

      if OutputCadence ~= 0 and n % OutputCadence == 0 then
         local fname = string.format('data/gravity-wave-%04d.h5', n)
         local outfile = hdf5.File(fname, 'w')
         outfile['prim' ] = P          [{{Ng,-Ng},nil}]
         outfile['exact'] = solution(t)[{{Ng,-Ng},nil}]
         outfile:close()
      end

      local start = os.clock()
      fish.grav1d_advance(dt)
      t = t + dt
      n = n + 1
      local kzps = 1e-3 * N / (os.clock() - start)

      if MessageCadence ~= 0 and n % MessageCadence == 0 then
         print(string.format("%05d: t=%3.2f dt=%2.1e %4.3fkz/s", n, t, dt, kzps))
      end
   end

   fish.grav1d_finalize()
   fish.del(scheme)
   fluids.descr_del(descr)

   local P0 = solution(t)[{{Ng,-Ng},nil}]:vector()
   local P1 = P          [{{Ng,-Ng},nil}]:vector()
   local L1 = 0.0

   for i=0,#P0-1,5 do
      L1 = L1 + math.abs(P1[i] - P0[i]) * dx
   end
   ErrorTable[N] = L1

   return P, solution(t)
end

local function resolution_study()
   for _,res in pairs{8,16,32,64,128,256,512,1024} do
      run(res)
   end
   util.pretty_print(ErrorTable)
   print("estimated convergence rate is",
	 math.log10(ErrorTable[1024]/ErrorTable[64]) / math.log10(1024/64))
   util.plot(ErrorTable, {'set logscale'})
end

local function runandplot()
   local P, Pexact = run(128)
   util.plot{['code' ]=P     [{{Ng,-Ng},{0,1}}]:table(),
	     ['exact']=Pexact[{{Ng,-Ng},{0,1}}]:table()}
end

runandplot()
