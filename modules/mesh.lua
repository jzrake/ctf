
local oo     = require 'class'
local array  = require 'array'
local fish   = require 'fish'
local flui   = require 'fluids'

local Block = oo.class('Block')
local FluidDescriptor = oo.class('FluidDescriptor')


function FluidDescriptor:__init__()
   local descr = flui.descr_new()
   flui.descr_setfluid(descr, flui.NRHYD)
   flui.descr_setgamma(descr, 1.4)
   flui.descr_seteos(descr, flui.EOS_GAMMALAW)
   self._descr = descr
end

function FluidDescriptor:__gc__()
   flui.descr_del(self._descr)
end

function FluidDescriptor:get_ncomp(which)
   return flui.descr_getncomp(self._descr, flui[which:upper()])
end

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
   local rank = #args.size
   local guard = args.guard or 0
   local parent = args.parent
   local descr = parent and parent._descr or FluidDescriptor()
   local totsize = { }

   fish.block_setdescr (block, descr._descr)
   fish.block_setrank  (block, rank)
   fish.block_setguard (block, guard)

   for i,N in ipairs(args.size) do
      fish.block_setsize  (block, i-1, N)
      fish.block_setrange (block, i-1, 0.0, 1.0)
      totsize[i] = N + 2 * guard
   end

   table.insert(totsize, descr:get_ncomp('primitive'))
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
      self._id                  = args.id
   else
      self._parent              = nil
      self._root                = self
      self._registry            = { }
      self._id                  = 0
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
   return string.format("<fish.block @ %s: [%s] level=%d id=%d>",
			string.sub(tostring(self._block), 11),
                        table.concat(self:size(), ' '), self:level(), self._id)
end
function Block:__index__(key)
   if type(key) == 'string' then
      return oo.getattrib(self, key)
   else
      return self:get_child_block(key)
   end
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
   if self._children[id] then
      error('child already exists')
   end
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

function Block:fill_guard()
   fish.block_fillguard(self._block)
end

function Block:fill_conserved()
   fish.block_fillconserved(self._block)
end

function Block:evolve(W, dt)
   fish.block_evolve(self._block, W:buffer(), dt)
end

function Block:time_derivative(scheme)
   fish.block_timederivative(self._block, scheme)
end

function Block:grid_spacing(opts)
   --
   -- Return the maximum wavespeed over the block
   --
   local opts = opts or { }
   local mode = opts.mode or 'local'
   if opts.mode == 'local' then
      return fish.block_gridspacing(self._block, 0)
   elseif opts.mode == 'smallest' then
      local dx_min = fish.block_gridspacing(self._block, 0)
      for b in self:walk() do
	 local dx = fish.block_gridspacing(b._block, 0)
	 if dx < dx_min then dx_min = dx end
      end
      return dx_min
   else
      error('mode must be one of [local, smallest]')
   end
end

function Block:total_states(opts)
-- **************************************************************************
   -- Return the total number of interior zones on this block
   --
   -- opts: {table} (optional)
   --
   --  + mode: string ('including_guard'), ['interior', 'including_guard']
   --  + recurse: boolean (false), accumulate size over descendant blocks
   --
   -- **************************************************************************
   --
   local opts = opts or { }
   local mode = (opts.mode or 'including_guard'):upper()
   local recurse = opts.recurse or false
   if opts.recurse then
      local ntot = 0
      for b in self:walk() do
	 ntot = ntot + fish.block_totalstates(self._block, fish[mode])
      end
      return ntot
   else
      return fish.block_totalstates(self._block, fish[mode])
   end
end

function Block:max_wavespeed()
   --
   -- Return the maximum wavespeed over the block and its descendants
   --
   local Amax = 0.0
   for b in self:walk() do
      local A = fish.block_maxwavespeed(self._block)
      if A > Amax then Amax = A end
   end
   return Amax
end

function Block:map(f)
   --
   -- Map the function f(x,y,z) over the coordinates of the block, filling in
   -- the primitive data values
   --
   local Nx, Ny, Nz = table.unpack(self:size())
   local Ng = self:guard()
   local Pvec = self._primitive:vector()
   for i=0,Nx+2*Ng-1 do
      local x = fish.block_positionatindex(self._block, 0, i)
      local Pi = f(x,0,0)
      Pvec[5*i + 0] = Pi[1]
      Pvec[5*i + 1] = Pi[2]
      Pvec[5*i + 2] = Pi[3]
      Pvec[5*i + 3] = Pi[4]
      Pvec[5*i + 4] = Pi[5]
   end
end

function Block:table(q, T, opts)
   --
   -- Return a table T[x] = P[q] for P in each interior zone over the block,
   -- where x is the coordinate (1d only)
   --
   local T = T or { }
   local Nx, Ny, Nz = table.unpack(self:size())
   local Ng = self:guard()
   local Pvec = self._primitive:vector()
   for i=0,Nx+2*Ng-1 do
      local x = fish.block_positionatindex(self._block, 0, i)
      T[x] = Pvec[5*i + q]
   end
   return T
end

function Block:next(id)
   if id == nil then
      return 0, self
   end

   for i=id,7 do
      if self._children[i] then
	 return 0, self._children[i]
      end
   end

   if self._parent then
      return self._parent:next(self._id+1)
   else
      return nil
   end
end

function Block:walk()
   local function next_block(s)
      s.id, s.b = s.b:next(s.id)
      if not s.first and s.b and s.b:level() == s.level then
	 return nil
      else
	 s.first = false
	 return s.b
      end
   end
   local state = {b=self,level=self:level(),first=true}
   return next_block, state
end

local function test1()
   local block0 = Block{ size={16}, guard=2 }
   local block1 = block0:add_child_block(0)
   local block2 = block1:add_child_block(0)
   local block3 = block1:add_child_block(1)

   block3:add_child_block(0):add_child_block(0)
   print(block2:get_neighbor_block(0, 'R'))

   block0:fill_guard()
   block2:fill_guard()
   block2:map(function(x) return {x+1,1,0,0,0} end)

   print(block2:total_states())
   print(block2:total_states{mode='interior'})
   print(block2:max_wavespeed())
end

local function test2()
   local mesh = Block{ size={16}, guard=2 }
   local N = 2
   for i=0,N-1 do
      mesh:add_child_block(i)
   end
   for i=0,N-1 do
      for j=0,N-1 do
	 mesh:get_child_block(i):add_child_block(j)
      end
   end
   for i=0,N-1 do
      for j=0,N-1 do
	 for k=0,N-1 do
	    local b = mesh:get_child_block(i):get_child_block(j):add_child_block(k)
	 end
      end
   end
   local n = 0
   for b in mesh:walk() do
      print(b)
      n = n + 1
   end
   print('there are '..n..' total blocks')
   print(mesh:total_states{recurse=true, mode='interior'} / n)

   local n = 0
   for b in mesh:get_child_block(0):walk() do
      print(b)
      n = n + 1
   end
   print('there are '..n..' total blocks')
   print(mesh:get_child_block(0):total_states{recurse=true, mode='interior'} / n)
   print(mesh:grid_spacing{mode='local'})
   print(mesh:grid_spacing{mode='smallest'})
end

if ... then -- if __name__ == "__main__"
   return {Block=Block}
else
   test1()
   test2()
   print(debug.getinfo(1).source, ": All tests passed")
end
