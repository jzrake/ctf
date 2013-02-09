
local oo     = require 'class'
local array  = require 'array'
local fish   = require 'fish'
local flui   = require 'fluids'

local Block = oo.class('Block')

function Block:__init__(args)
   local block = fish.block_new()
   local descr = flui.descr_new()

   flui.descr_setfluid(descr, flui.NRHYD)
   flui.descr_setgamma(descr, 1.4)
   flui.descr_seteos(descr, flui.EOS_GAMMALAW)

   fish.block_setdescr (block, descr)
   fish.block_setrank  (block, #args.size)
   fish.block_setguard (block, args.guard or 0)

   for i,N in ipairs(args.size) do
      fish.block_setsize  (block, i-1, N)
      fish.block_setrange (block, i-1, 0.0, 1.0)
   end

   fish.block_allocate(block)

   self._children = { }
   self._rank = #args.size
   self._block = block
   self._descr = descr
end

function Block:__gc__(args)
   fish.block_del(self._block)
   flui.descr_del(self._descr)
end

function Block:__tostring__(args)
   return string.format("<fish.block: [%s] level %d>",
			table.concat(self:size(), ' '), self:level())
end

function Block:level()
   return self._parent and self._parent:level() + 1 or 0
end

function Block:size()
   local S = { }
   for i=1,self._rank do
      S[i] = fish.block_getsize(self._block, i-1)
   end
   return S
end

function Block:guard()
   return fish.block_getguard(self._block)
end

function Block:create_child_block(id)
   local child = Block{ size=self:size(),
			guard=self:guard() }
   fish.block_setchild(self._block, id, child._block)
   self._children[id] = child
   child._parent = self
   return child
end


local function test1()
   local block0 = Block{ size={16,16,16}, guard=2 }
   local block1 = block0:create_child_block(0)
   local block2 = block1:create_child_block(0)
   print(block2)
end


if ... then -- if __name__ == "__main__"
   return {Block=Block}
else
   test1()
   print(debug.getinfo(1).source, ": All tests passed")
end
