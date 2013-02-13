
local fish   = require 'fish'
local fluids = require 'fluids'
local array  = require 'array'

local function gzone()
   local Ng = 5
   for i=0,Ng-1 do
      local n0 = math.floor((i + Ng - 1)/2)
      --print(i, n0, n0+1)
   end

   for i=0,Ng-1 do
      local n0 = -1 + math.floor((i+1)/2)
      print(i, n0, n0+1)
   end
end

local function testinterp()
   local Nx = 64
   local Ng = 4
   local B0 = fish.block_new()
   local BL = fish.block_new()
   local BR = fish.block_new()
   local bufs = { }
   local descr = fluids.descr_new()

   fluids.descr_setfluid(descr, fluids.NRHYD)
   fluids.descr_setgamma(descr, 1.4)
   fluids.descr_seteos(descr, fluids.EOS_GAMMALAW)

   for _,B in pairs{B0, BL, BR} do

      local P = array.array{Nx + 2*Ng, 5}
      local Pvec = P:vector()

      fish.block_setdescr(B, descr)
      fish.block_setguard(B, Ng)
      fish.block_setrank(B, 1)
      fish.block_setsize(B, 0, Nx)
      fish.block_allocate(B)
      fish.block_mapbuffer(B, P:buffer(), fluids.PRIMITIVE)

      bufs[B] = P
   end

   fish.block_setrange(B0, 0, 0.0, 1.0)

   fish.block_setchild(B0, 0, BL) -- id=0
   fish.block_setchild(B0, 1, BR) -- id=1
   --fish.block_setneighbor(BL, 0, fish.RIGHT, BR)
   --fish.block_setneighbor(BR, 0, fish.LEFT, BL)

   for _,B in pairs{B0, BL, BR} do

      local Pvec = bufs[B]:vector()
      for i=0,Nx+2*Ng-1 do
	 local x = fish.block_positionatindex(B, 0, i)
	 Pvec[i*5] = math.sin(2 * math.pi * x)
      end

   end

   fish.block_fillguard(BL)
   fish.block_fillguard(BR)

   local names = {[B0]='B0.dat', [BL]='BL.dat', [BR]='BR.dat'}

   for _,B in pairs{B0, BL, BR} do
      local f = io.open(names[B], 'w')
      local Pvec = bufs[B]:vector()
      for i=0,Nx+2*Ng-1 do
	 local x = fish.block_positionatindex(B, 0, i)
	 f:write(string.format('%f %f\n', x, Pvec[5*i]))
      end
      f:close()
   end

   fluids.descr_del(descr)
   fish.block_del(B0)
   fish.block_del(BL)
   fish.block_del(BR)
end

--gzone()
testinterp()
