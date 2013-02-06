
local array  = require 'array'
local fish   = require 'fish'
local fluids = require 'fluids'

local function test1()
   local B = fish.block_new()
   local err = fish.block_setrank(B, 4)
   assert(err < 0)
   assert(fish.block_geterror(B) == "rank must be 1, 2, or 3")
   fish.block_setrank(B, 3)
   assert(not fish.block_geterror(B))
   fish.block_del(B)
end

local function test2()
   local D = fluids.descr_new()
   local B = fish.block_new()

   fluids.descr_setfluid(D, fluids.NRHYD)
   fluids.descr_setgamma(D, 1.4)
   fluids.descr_seteos(D, fluids.EOS_GAMMALAW)

   fish.block_setrank(B, 2)
   fish.block_setsize(B, 0, 128)
   fish.block_setsize(B, 1, 256)
   assert(fish.block_getrank(B) == 2)
   assert(fish.block_getsize(B, 0) == 128)
   assert(fish.block_getsize(B, 1) == 256)
   fish.block_getsize(B, 2)
   assert(fish.block_geterror(B) ==
	  "argument 'dim' must be smaller than the rank of the block")
   fish.block_allocate(B)
   assert(fish.block_geterror(B) ==
	  "block's fluid descriptor must be set before allocating")
   fish.block_setdescr(B, D)
   assert(fish.block_allocate(B) == 0)
   assert(fish.block_allocate(B) == 0) -- reallocates the fluid states
   fish.block_del(B)
   fluids.descr_del(D)
end

local function test3()
   local x0 = array.vector(1)
   local x1 = array.vector(1)
   local B = fish.block_new()
   fish.block_setrank(B, 1)
   fish.block_setsize(B, 0, 20)
   fish.block_setrange(B, 0, -1.0, 1.0)
   fish.block_getrange(B, 0, x0:pointer(), x1:pointer())
   fish.block_allocate(B)
   assert(x0[0] == -1.0)
   assert(x1[0] ==  1.0)
   assert(fish.block_gridspacing(B, 0) == 0.1)
   fish.block_del(B)
end

local function test4()
   local B0 = fish.block_new()
   local BL = fish.block_new()
   local BR = fish.block_new()

   for _,B in pairs{B0, BL, BR} do
      fish.block_setrank(B, 3)
      fish.block_setsize(B, 0, 16)
      fish.block_setsize(B, 1, 16)
      fish.block_setsize(B, 2, 16)
      fish.block_setrange(B, 0, -1.0, 1.0)
      fish.block_setrange(B, 1, -1.0, 1.0)
      fish.block_setrange(B, 2, -1.0, 1.0)
      fish.block_allocate(B)
   end

   local BchildL = fish.block_light()
   local BchildR = fish.block_light()
   local x0 = array.vector(1)
   local x1 = array.vector(1)

   fish.block_setchild(B0, 0, BL) -- id=0
   fish.block_setchild(B0, 3, BR) -- id=3
   fish.block_getchild(B0, 0, BchildL) -- id=0
   fish.block_getchild(B0, 3, BchildR) -- id=3

   assert(fish.block_light(BchildL) == fish.block_light(BL))
   assert(fish.block_light(BchildR) == fish.block_light(BR))

   fish.block_getrange(BchildL, 0, x0:pointer(), x1:pointer())
   assert(x0[0] == -1.0)
   assert(x1[0] ==  0.0)
   fish.block_getrange(BchildL, 1, x0:pointer(), x1:pointer())
   assert(x0[0] == -1.0)
   assert(x1[0] ==  0.0)
   fish.block_getrange(BchildL, 2, x0:pointer(), x1:pointer())
   assert(x0[0] == -1.0)
   assert(x1[0] ==  0.0)

   fish.block_getrange(BchildR, 0, x0:pointer(), x1:pointer())
   assert(x0[0] ==  0.0)
   assert(x1[0] ==  1.0)
   fish.block_getrange(BchildR, 1, x0:pointer(), x1:pointer())
   assert(x0[0] ==  0.0)
   assert(x1[0] ==  1.0)
   fish.block_getrange(BchildR, 2, x0:pointer(), x1:pointer())
   assert(x0[0] == -1.0)
   assert(x1[0] ==  0.0)

   fish.block_del(B0)
   fish.block_del(BL)
   fish.block_del(BR)
end

test1()
test2()
test3()
test4()
