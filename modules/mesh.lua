
local oo     = require 'class'
local array  = require 'array'
local fish   = require 'fish'
local flui   = require 'fluids'

local Block = oo.class('Block')

function Block:__init__(args)
   -- **************************************************************************
   -- args: {table}
   --
   --  + size: table, e.g. {16,32,32}
   --  + guard: number, of guard zones
   --  + parent: block, if nil create new root block
   --  + id: number, needed if parent is given
   --
   -- **************************************************************************
   local block = fish.block_new()
   local descr = flui.descr_new()
   local rank = #args.size
   local guard = args.guard or 0
   local parent = args.parent
   local totsize = { }

   flui.descr_setfluid(descr, flui.NRHYD)
   flui.descr_setgamma(descr, 1.4)
   flui.descr_seteos(descr, flui.EOS_GAMMALAW)

   fish.block_setdescr (block, descr)
   fish.block_setrank  (block, rank)
   fish.block_setguard (block, guard)

   for i,N in ipairs(args.size) do
      fish.block_setsize  (block, i-1, N)
      fish.block_setrange (block, i-1, 0.0, 1.0)
      totsize[i] = N + 2 * guard
   end

   table.insert(totsize, flui.descr_getncomp(descr, flui.PRIMITIVE))
   local primitive = array.array(totsize)

   fish.block_allocate  (block)
   fish.block_mapbuffer (block, primitive:buffer(), flui.PRIMITIVE)

   if parent then
      if not args.id then
	 error("id must be given if creating a child block")
      end
      fish.block_setchild(parent._block, args.id, block)
      parent._children[args.id] = self

      self._parent              = parent
      self._root                = parent._root
      self._registry            = parent._registry
   else
      self._parent              = nil
      self._root                = self
      self._registry            = { }
   end

   self._primitive = primitive
   self._block = block
   self._descr = descr
   self._children = { }
   self._rank = rank
   self._registry[fish.block_light(block)] = self
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
   --
   -- Return the level depth of the block below root
   --
   return self._parent and self._parent:level() + 1 or 0
end

function Block:size()
   --
   -- Return a numeric table containing the block size (not including guard
   -- zones)
   --
   local S = { }
   for i=1,self._rank do
      S[i] = fish.block_getsize(self._block, i-1)
   end
   return S
end

function Block:guard()
   --
   -- Return the number of guard zones
   --
   return fish.block_getguard(self._block)
end

function Block:add_child_block(id)
   --
   -- Create and return a new child block with index `id`
   --
   local child = Block{ size=self:size(),
			guard=self:guard(),
			parent=self,
			id=id }
   return child
end

function Block:get_child_block(id)
   --
   -- Return the child block at index `id`
   --
   return self._children[id]
end

function Block:get_neighbor_block(dim, dir)
   -- **************************************************************************
   -- args:
   --
   --  + dim: number, 0, 1, or 2
   --  + dir: 'L' or 'R' (left or right neighbor)
   --
   -- **************************************************************************
   local d = ({L=fish.LEFT, R=fish.RIGHT})[dir]
   local c = fish.block_getneighbor(self._block, dim, d)
   return self._registry[c]
end


local function test1()
   local block0 = Block{ size={16,16,16}, guard=2 }
   local block1 = block0:add_child_block(0)
   local block2 = block1:add_child_block(0)
   local block3 = block1:add_child_block(1)

   print(block2)
   print(block2:get_neighbor_block(0, 'R'))

   for k,v in pairs(block0._registry) do
      print(k,v)
   end
end


if ... then -- if __name__ == "__main__"
   return {Block=Block}
else
   test1()
   print(debug.getinfo(1).source, ": All tests passed")
end
