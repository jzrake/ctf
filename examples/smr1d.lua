local fish = require 'fish'

local function test1()
   local B = fish.block_new()
   local err = fish.block_setrank(B, 4)
   assert(err == fish.ERROR_BADARG)
   assert(fish.block_geterror(B) == "rank must be 1, 2, or 3")
   fish.block_del(B)
end

local function test2()
   local B = fish.block_new()
   fish.block_setrank(B, 2)
   fish.block_setsize(B, 0, 128)
   fish.block_setsize(B, 1, 256)
   assert(fish.block_getrank(B) == 2)
   assert(fish.block_getsize(B, 0) == 128)
   assert(fish.block_getsize(B, 1) == 256)
   fish.block_getsize(B, 2)
   assert(fish.block_geterror(B) ==
	  "argument 'dim' must be smaller than the rank of the block")
   fish.block_del(B)
end

test1()
test2()
