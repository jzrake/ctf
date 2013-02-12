
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

function FluidDescriptor:fluid()
   local E = { }
   for k,v in pairs(flui) do
      if type(v) == 'number' then
	 E[v] = k
      end
   end
   local F = array.vector(1, 'int')
   flui.descr_getfluid(self._descr, F:buffer())
   return E[F[0]]:lower()
end

function Block:__init__(args)
   -- **************************************************************************
   -- args: {table}
   --
   --  + size: table, e.g. {16,32,32}
   --  + guard: number, of guard zones
   --  + parent: block, if nil create new root block
   --  + id: number, needed if parent is given
   --  + dummy: boolean, false by default: block will not be allocated
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

   if not args.dummy then
      local primitive = array.array(totsize)
      fish.block_allocate  (block)
      fish.block_mapbuffer (block, primitive:buffer(), flui.PRIMITIVE)
      self._primitive = primitive
   else
      self._primitive = { }
   end

   self._block = block
   self._descr = descr
   self._children = { }
   self._boundaryL = { }
   self._boundaryR = { }
   self._rank = rank

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

      self:set_boundary_block(0, 'L', self) -- periodic BC's on root
      self:set_boundary_block(0, 'R', self)
   end

   self._registry[fish.block_light(block)] = self
end

function Block:__gc__(args)
   fish.block_del(self._block)
end

function Block:__tostring__(args)
   return string.format("<fish.block @ %s: [%s] level=%d id=%d>",
			string.sub(tostring(self._block), 11),
                        table.concat(self:size(), ' '), self:level(), self._id)
end

function Block:__index__(key)
   if type(key) == 'string' then
      return ( -- These are the class read-only attribute
	 {descr=self._descr,
	  primitive=self._primitive})[key] or oo.getattrib(self, key)
   else
      return self:child_block(key)
   end
end

function Block:level()
   --
   -- Return the level depth of the block below root
   --
   return fish.block_level(self._block)
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

function Block:dummy()
   return fish.block_allocated(self._block) == 0 and true or false
end

function Block:add_child_block(id, opts)
   --
   -- Create and return a new child block with index `id`
   --
   if self._children[id] then
      error('child already exists')
   end
   local child = Block{ size=self:size(),
                        guard=self:guard(),
                        parent=self,
                        id=id,
			dummy=(opts or { }).dummy }
   return child
end

function Block:child_block(id)
   --
   -- Return the child block at index `id`
   --
   return self._children[id]
end

function Block:has_children()
   return next(self._children) and true or false
end

function Block:neighbor_block(dim, dir)
   -- **************************************************************************
   -- args:
   --
   --  + dim: number, 0, 1, or 2
   --  + dir: 'L' or 'R' (left or right neighbor)
   --
   -- **************************************************************************
   local d = ({L=fish.LEFT, R=fish.RIGHT})[dir]
   local c = fish.block_neighbor(self._block, dim, d)
   return self._registry[c]
end

function Block:get_boundary_block(dim, dir)
   --
   -- Set a boundary block to over-ride mesh connectivity. Arguments are the
   -- same as Block:neighbor_block.
   --
   local d = ({L=fish.LEFT, R=fish.RIGHT})[dir]
   local c = fish.block_getboundaryblock(self._block, dim, d)
   return self._registry[c]
end

function Block:set_boundary_block(dim, dir, block)
   --
   -- Get the boundary block. Arguments are the same as Block:neighbor_block,
   -- but with the boundary block given as the third argument.
   --
   local d = ({L=fish.LEFT, R=fish.RIGHT})[dir]
   self['_boundary'..dir][dim] = block
   fish.block_setboundaryblock(self._block, dim, d, block._block)
end

function Block:fill_conserved()
   fish.block_fillconserved(self._block)
end

function Block:fill_guard()
   fish.block_fillguard(self._block)
end

function Block:evolve(W, dt)
   fish.block_evolve(self._block, W:buffer(), dt)
end

function Block:time_derivative(scheme)
   fish.block_timederivative(self._block, scheme)
end

function Block:project()
   fish.block_project(self._block)
end

function Block:fill()
   for _,c in pairs(self._children) do
      if c:has_children() then
	 c:fill()
      end
      if not self:dummy() then
	 c:project()
      end
   end
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
      local A = fish.block_maxwavespeed(b._block)
      if A > Amax then Amax = A end
   end
   return Amax
end

function Block:map(f, opts)
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
   for i=Ng,Nx+Ng-1 do
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

function Block:walk(opts)
   --
   -- Return an iterator to traverse the tree of blocks depth-first. Unless
   -- opts.include_dummy == true, inactive (dummy) blocks are skipped.
   --
   local opts = opts or { }
   local function next_block(s)
      repeat -- skip dummy blocks in walk
	 s.id, s.b = s.b:next(s.id)
      until opts.include_dummy or not (s.b and s.b:dummy())
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
   assert(block2:neighbor_block(0, 'R') == block3)

   block0:fill_guard()
   block2:fill_guard()
   block2:map(function(x) return {x+1,1,0,0,0} end)

   assert(block2:total_states() == 20)
   assert(block2:total_states{mode='interior'} == 16)
   assert(math.abs(block2:max_wavespeed() - 1.1973303637677) < 1e-10)
end

local function test2()
   local mesh = Block{ size={16}, guard=2 }
   local N = 2
   for i=0,N-1 do
      mesh:add_child_block(i)
   end
   for i=0,N-1 do
      for j=0,N-1 do
	 mesh[i]:add_child_block(j)
      end
   end
   for i=0,N-1 do
      for j=0,N-1 do
	 for k=0,N-1 do
	    local b = mesh[i][j]:add_child_block(k)
	 end
      end
   end
   local n = 0
   for b in mesh:walk() do n = n + 1 end
   assert(n == 15)
   assert(mesh:total_states{recurse=true, mode='interior'} / n == 16)

   local n = 0
   for b in mesh[0]:walk() do n = n + 1 end

   assert(n == 7)
   assert(mesh:grid_spacing{mode='local'} == 0.0625)
   assert(mesh:grid_spacing{mode='smallest'} == 0.0078125)
end

local function test3()
   local mesh = Block { size={32} }
   assert(mesh.descr:fluid() == 'nrhyd')

   mesh:set_boundary_block(0, 'L', mesh)
   mesh:set_boundary_block(0, 'R', mesh)
   assert(mesh:get_boundary_block(0, 'L') == mesh)
   assert(mesh:get_boundary_block(0, 'R') == mesh)

   mesh:add_child_block(0)
   local err, msg = pcall(mesh.set_boundary_block, mesh, 0, 'L', mesh[0])
   assert(not err)

   local dummy = Block { size={32,32}, dummy=true }
   dummy      :add_child_block(0, {dummy=true})
   dummy[0]   :add_child_block(0, {dummy=false})
   dummy[0][0]:add_child_block(0, {dummy=false})
   assert      (dummy:dummy())
   assert(   dummy[0]:dummy())
   assert(not mesh[0]:dummy())

   local F = function(x) return {x+1,1,0,0,0} end
   for b in dummy:walk() do
      b:map(F)
   end
   dummy:fill()
end

local function test4()
   local mesh = Block{ size={64}, guard=3 }
   local block1 = mesh:add_child_block(0)
   local block2 = mesh:add_child_block(1)

   local F = function(x) return {x+1,1,0,0,0} end
   block1:map(F)
   block2:map(F)
   assert(mesh:has_children())
   assert(not block1:has_children())

   mesh:fill()
   --block1:project()
   --block2:project()

   local P = mesh.primitive[{{3,-3},{0,1}}]:vector()
   for i=0,#P-1 do
      local x = fish.block_positionatindex(mesh._block, 0, i+3)
      assert(P[i] == F(x)[1])
   end
end

if ... then -- if __name__ == "__main__"
   return {Block=Block}
else
   test1()
   test2()
   test3()
   test4()
   print(debug.getinfo(1).source, ": All tests passed")
end
