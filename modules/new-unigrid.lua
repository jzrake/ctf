
local MPI    = require 'MPI'
local array  = require 'array'
local oo     = require 'class'
local cow    = require 'cow'
local hdf5   = require 'lua-hdf5.LuaHDF5'
local util   = require 'util'

local UnigridDomain    = oo.class('UnigridDomain')
local UnigridDataField = oo.class('UnigridDataField')


--------------------------------------------------------------------------------
-- UnigridDomain implementation
--------------------------------------------------------------------------------
function UnigridDomain:__init__(size, guard)
   self._domain = cow.domain_new()
   self._comm = next(MPI) and MPI.Comm() or nil
   cow.domain_setndim(self._domain, #size)
   cow.domain_setsize(self._domain, 0, size[1] or 1)
   cow.domain_setsize(self._domain, 1, size[2] or 1)
   cow.domain_setsize(self._domain, 2, size[3] or 1)
   cow.domain_setguard(self._domain, guard or 0)
   cow.domain_commit(self._domain)
   cow.domain_getcomm(self._domain, self._comm)
end

function UnigridDomain:__gc__(size, guard)
   cow.domain_del(self._domain)
end

function UnigridDomain:new_data_field(members)
   return UnigridDataField(self, members)
end

function UnigridDomain:size(dim, which)
   --
   -- dim: 0, 1, 2, 'all' (default), or 'shape'
   -- which: 'global' (default), 'local-interior' or 'local-including-guard'
   --
   local dim = dim or 'all'
   local which = which or 'global'
   local d, func

   if dim == 'all' then
      d = cow.ALL_DIMS
   elseif dim == 'shape' then
      d = 'shape'
   elseif dim == 0 or dim == 1 or dim == 2 then
      d = dim
   else
      error("dim: 0, 1, 2, 'all', or 'shape'")
   end

   if which == 'global' then
      func = cow.domain_getnumglobalzones
   elseif which == 'interior' then
      func = cow.domain_getnumlocalzonesinterior
   elseif which == 'local-including-guard' then
      func = cow.domain_getnumlocalzonesincguard
   else
      error("which: 'global', 'local-interior' or 'local-including-guard'")
   end

   if type(d) == 'number' then
      return func(self._domain, d)
   else
      local ret = { }
      for i=1,cow.domain_getndim(self._domain) do
	 ret[i] = func(self._domain, i-1)
      end
      return ret
   end
end

function UnigridDomain:get_communicator(members)
   return self._comm
end

--------------------------------------------------------------------------------
-- UnigridDataField implementation
--------------------------------------------------------------------------------
function UnigridDataField:__init__(domain, members)
   self._dfield = cow.dfield_new()
   cow.dfield_setdomain(self._dfield, domain._domain)
   for _,member in pairs(members) do
      cow.dfield_addmember(self._dfield, member)
   end
   local size = domain:size('shape', 'local-including-guard')

   size[#size + 1] = #members
   self._buffer = array.array(size)

   cow.dfield_setdatabuffer(self._dfield, self._buffer:buffer())
   cow.dfield_commit(self._dfield)
end

function UnigridDataField:__gc__(size, guard)
   cow.dfield_del(self._dfield)
end

function UnigridDataField:read(fname, gname)
   cow.dfield_setname(self._dfield, gname)
   cow.dfield_read(self._dfield, fname)
end

function UnigridDataField:write(fname, gname)
   cow.dfield_setname(self._dfield, gname)
   cow.dfield_write(self._dfield, fname)
end

return { UnigridDomain=UnigridDomain,
	 UnigridDataField=UnigridDataField }
