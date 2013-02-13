
--------------------------------------------------------------------------------
--
--               High level Lua wrappers for HDF5 library
--
--------------------------------------------------------------------------------
local hdf5 = { } -- module table
--------------------------------------------------------------------------------

local H5      = require 'HDF5'
local buffer  = require 'buffer'
local array   = require 'array'
local oo      = require 'class'
local hp0     = H5.H5P_DEFAULT

if not next(H5) then
   print('warning! HDF5 was imported but is not available')
   return { }
end

local Base      = oo.class('HDF5_BaseClass')
local Indexable = oo.class('HDF5_BaseClass', Base)

hdf5.File      = oo.class('File', Indexable)
hdf5.Group     = oo.class('Group', Indexable)
hdf5.DataSet   = oo.class('DataSet', Base)
hdf5.DataSpace = oo.class('DataSpace', Base)
hdf5.DataType  = oo.class('DataType', Base)


local function mpio_stats(dxpl, mpio)
   -----------------------------------------------------------------------------
   -- For HDF5/MPIO only. Call will simply have no result without parallel HDF5
   -- support. On HDF5 versions < 1.8.10, the [global/local]_no_collective_cause
   -- variables will not be evaluated. On HDF5 versions less than 1.8.8 none of
   -- these functions work.
   -----------------------------------------------------------------------------
   if not H5.H5_VERSION_GE(1,8,8) then return end
   local Lcause, Gcause = H5.H5Pget_mpio_no_collective_cause(dxpl)
   local L = { }
   local G = { }
   if Lcause and Gcause then
      if Lcause == H5.H5D_MPIO_COLLECTIVE then L[1] = 'COLLECTIVE' end
      if Gcause == H5.H5D_MPIO_COLLECTIVE then G[1] = 'COLLECTIVE' end
      local s = 'H5D_MPIO_'
      for _,c in ipairs({'SET_INDEPENDENT',
			 'DATATYPE_CONVERSION',
			 'DATA_TRANSFORMS',
			 'SET_MPIPOSIX',
			 'NOT_SIMPLE_OR_SCALAR_DATASPACES',
			 'POINT_SELECTIONS',
			 'NOT_CONTIGUOUS_OR_CHUNKED_DATASET',
			 'FILTERS'}) do
	 if bit32.band(Lcause, H5[s..c]) ~= 0 then table.insert(L, c) end
	 if bit32.band(Gcause, H5[s..c]) ~= 0 then table.insert(G, c) end
      end
   end
   mpio.local_no_collective_cause = table.concat(L, ', ')
   mpio.global_no_collective_cause = table.concat(G, ', ')
   mpio.actual_chunk_opt_mode = (
      {[H5.H5D_MPIO_NO_CHUNK_OPTIMIZATION]='NO_CHUNK_OPTIMIZATION',
       [H5.H5D_MPIO_MULTI_CHUNK]='MULTI_CHUNK',
       [H5.H5D_MPIO_MULTI_CHUNK_NO_OPT]="MULTI_CHUNK_NO_OPT:",
       [H5.H5D_MPIO_LINK_CHUNK]='LINK_CHUNK'}
   )[H5.H5Pget_mpio_actual_chunk_opt_mode(dxpl)]
   mpio.actual_io_mode = (
      {[H5.H5D_MPIO_NO_COLLECTIVE]="NO_COLLECTIVE",
       [H5.H5D_MPIO_CHUNK_INDEPENDENT]="CHUNK_INDEPENDENT",
       [H5.H5D_MPIO_CHUNK_COLLECTIVE]="CHUNK_COLLECTIVE",
       [H5.H5D_MPIO_CHUNK_MIXED]="CHUNK_MIXED",
       [H5.H5D_MPIO_CONTIGUOUS_COLLECTIVE]="CONTIGUOUS_COLLECTIVE"}
   )[H5.H5Pget_mpio_actual_io_mode(dxpl)]
end


--------------------------------------------------------------------------------
-- Base classes for meta-table and methods wrapping hid_t objects
--------------------------------------------------------------------------------
function Base:__index__()
   if self._hid == 0 then -- object is closed
      return nil
   end
end

function Base:__newindex__()
   if self._hid == 0 then -- object is closed
      error("cannot assign to closed object")
   end
end

function Base:__gc__()
   self:close()
end

function Base:close()
   if self._open_objects then
      for k,v in pairs(self._open_objects) do
	 v:close()
	 self._open_objects[k] = nil
      end
   end
   if self._parent then
      self._parent._open_objects[self._name] = nil
   end
   if self._hid ~= 0 then
      self._close(self._hid)
      self._hid = 0
   else
      -- object is already closed
   end
end

function Base:class()
   return self._type
end

function Base:name()
   return self._name
end


--------------------------------------------------------------------------------
-- Indexable methods
--------------------------------------------------------------------------------
function Indexable:path(key)
   if not self._parent then
      return '/' .. self._name
   else
      return self._parent:path() .. '/' .. self._name
   end
end

function Indexable:keys()
   local link_names = { }
   local function f(name)
      table.insert(link_names, name)
   end
   local idx = H5.new_hsize_t_arr{0}
   H5.H5Literate(self._hid, H5.H5_INDEX_NAME, H5.H5_ITER_NATIVE, idx, f)
   return link_names
end

function Indexable:require_group(name)
   return hdf5.Group(self, name)
end

function Indexable:__index__(key)
   if self._hid == 0 then -- object is closed
      return nil
   end
   if rawget(self, '__dict__')._open_objects[key] then
      return rawget(self, '__dict__')._open_objects[key]
   end
   if not H5.H5Lexists(self._hid, key, hp0) then
      return nil
   end
   local info = H5.new_H5O_info_t()
   H5.H5Oget_info_by_name(self._hid, key, info, hp0)
   if info.type == H5.H5O_TYPE_GROUP then
      return hdf5.Group(self, key)
   elseif info.type == H5.H5O_TYPE_DATASET then
      return hdf5.DataSet(self, key)
   else
      error("Indexable:object %s/%s has an unsupported type", self._name, key)
   end
end

function Indexable:__newindex__(key, value)
   if self._hid == 0 then -- object is closed
      error("Indexable:cannot assign to closed object")
   end
   if type(value) == 'string' then
      local targ = self._hid
      local data = buffer.new_buffer(value)
      local fspc = H5.H5Screate(H5.H5S_SCALAR)
      local strn = H5.H5Tcopy(H5.H5T_C_S1)

      H5.H5Tset_size(strn, #data)
      local dset = H5.H5Dcreate2(targ, key, strn, fspc, hp0, hp0, hp0)

      H5.H5Dwrite(dset, strn, fspc, fspc, hp0, data)
      H5.H5Dclose(dset)
      H5.H5Tclose(strn)
      H5.H5Sclose(fspc)

   elseif type(value) == 'number' then
      local targ = self._hid
      local fspc = H5.H5Screate(H5.H5S_SCALAR)
      local data = array.vector{value}
      local dset = H5.H5Dcreate2(targ, key, H5.H5T_NATIVE_DOUBLE, fspc,
				 hp0, hp0, hp0)
      H5.H5Dwrite(dset, H5.H5T_NATIVE_DOUBLE, fspc, fspc, hp0, data:buffer())
      H5.H5Dclose(dset)
      H5.H5Sclose(fspc)

   elseif value.buffer and value.dtype and value.selection then
      --------------------------------------------------------------------------
      -- If the buffer, dtype, and selection methods are given, assume their
      -- behavior is like that of array.view, and we can write to a buffer
      -- automatically.
      --------------------------------------------------------------------------
      local start, stride, count, block = value:selection()
      local mspace = hdf5.DataSpace(value:extent())
      mspace:select_hyperslab(start, stride, count, block)

      local dset = hdf5.DataSet(self, key, 'w',
				{dtype=value:dtype(), shape=count})
      dset:write(value:buffer(), mspace, nil)
   else
      error("DataSet:unrecognized type for writing")
   end
end
function Indexable:__pairs__()
   local p = { }
   for k,v in pairs(self:keys()) do
      p[v] = self[v]
   end
   return pairs(p)
end


--------------------------------------------------------------------------------
-- HDF5 File class methods
--------------------------------------------------------------------------------
function hdf5.File:__tostring__()
   if self._hid ~= 0 then
      return string.format("<HDF5 file: \"%s\" (mode %s)>", self._name,
			   self._mode)
   else
      return "<Closed HDF5 file>"
   end
end


--------------------------------------------------------------------------------
-- HDF5 Group class methods
--------------------------------------------------------------------------------
function hdf5.Group:__tostring__()
   if self._hid ~= 0 then
      return string.format("<HDF5 group: \"%s\">", self._name)
   else
      return "<Closed HDF5 group>"
   end
end


--------------------------------------------------------------------------------
-- HDF5 DataSet class methods
--------------------------------------------------------------------------------
function hdf5.DataSet:read(buf, mspace, fspace)
   -----------------------------------------------------------------------------
   -- Read from the data set into `buf` according to source `fspace` and
   -- destination `mspace`, both of which default to all of the data set
   -- extent. If `buf` is also absent then it is created with the smallest
   -- possible size and returned.
   -----------------------------------------------------------------------------
   local fspace = fspace or self:get_space()
   local mspace = mspace or self:get_space()
   local htype = self:get_type()
   local bytes = mspace:get_select_npoints() * htype:get_size() -- source size
   local buf = buf or buffer.new_buffer(bytes)
   if bytes > #buf then
      error("data space selection is too large for buffer")
   end
   local dxpl = H5.H5Pcreate(H5.H5P_DATASET_XFER)
   local mpio_mode = H5['H5FD_MPIO_'..(self._mpio.requested_mode or '')]
   if mpio_mode and H5.H5Pset_dxpl_mpio then
      H5.H5Pset_dxpl_mpio(dxpl, mpio_mode)
   end
   local err = H5.H5Dread(self._hid, htype._hid, mspace._hid, fspace._hid, dxpl,
			  buf)
   if mpio_mode then
      mpio_stats(dxpl, self._mpio)
   end
   H5.H5Pclose(dxpl)
   if #err < 0 then error("DataSet:read") end
   return buf
end

function hdf5.DataSet:write(buf, mspace, fspace)
   -----------------------------------------------------------------------------
   -- Write from the data set into `buf` according to source `fspace` and
   -- destination `mspace`, both of which default to all of the data set
   -- extent.
   -----------------------------------------------------------------------------
   local fspace = fspace or self:get_space()
   local mspace = mspace or self:get_space()
   local htype = self:get_type()
   local bytes = fspace:get_select_npoints() * htype:get_size() -- dest size
   if bytes > #buf then
      error("DataSet:data space selection is too small for buffer")
   end
   local dxpl = H5.H5Pcreate(H5.H5P_DATASET_XFER)
   local mpio_mode = H5['H5FD_MPIO_'..(self._mpio.requested_mode or '')]
   if mpio_mode and H5.H5Pset_dxpl_mpio then
      H5.H5Pset_dxpl_mpio(dxpl, mpio_mode)
   end
   local err = H5.H5Dwrite(self._hid, htype._hid, mspace._hid, fspace._hid, dxpl,
			   buf)
   if mpio_mode then
      mpio_stats(dxpl, self._mpio)
   end
   H5.H5Pclose(dxpl)
   if #err < 0 then error("DataSet:write") end
end

function hdf5.DataSet:value()
   -----------------------------------------------------------------------------
   -- Read internal data into an object whose Lua type is inferred from the HDF5
   -- type. Might throw an error if it can't find an appropriate Lua type.
   -----------------------------------------------------------------------------
   local space = self:get_space()
   local start = { }
   local tstr = self:get_type():type_string()
   local tcls = self:get_type():type_class()
   if tcls == 'string' then return tostring(self:read())
   elseif tcls == 'float' then
      local extent = space:get_extent()
      local dtype = buffer[self:get_type():type_string()]
      if #extent == 0 then -- scalar data set
	 return buffer.get_typed(self:read(), dtype, 0)
      else
	 return array.view(self:read(), tstr, extent)
      end
   else error("DataSet:could not infer a Lua type from the data set")
   end
end

function hdf5.DataSet:get_space()
   return hdf5.DataSpace(self)
end

function hdf5.DataSet:get_chunk()
   local dcpl = H5.H5Dget_create_plist(self._hid)
   if H5.H5Pget_layout(dcpl) ~= H5.H5D_CHUNKED then
      return { }
   end
   local rank = #self:get_space():get_extent()
   local lchunk = { }
   for i=1,rank do lchunk[i] = 0 end
   local hchunk = H5.new_hsize_t_arr(lchunk)
   local ret = H5.H5Pget_chunk(dcpl, rank, hchunk)
   if ret < 0 then error('DataSet:get_chunk failed') end
   H5.H5Pclose(dcpl)
   for i=1,rank do lchunk[i] = hchunk[i-1] end
   return lchunk
end

function hdf5.DataSet:set_extent(extent)
   if not self:get_chunk() then
      error("DataSet:must be chunked to change its extent")
   end
   local rank = #self:get_space():get_extent()
   if rank ~= #extent then
      error("DataSet:new extent must match data set rank")
   end
   local hextent = H5.new_hsize_t_arr(extent)
   local err = H5.H5Dset_extent(self._hid, hextent)
   if #err < 0 then
      error("DataSet:set_extent failed")
   end
end

function hdf5.DataSet:get_type()
   -----------------------------------------------------------------------------
   -- Return a DataType class representing the DataSet's HDF5 type.
   -----------------------------------------------------------------------------
   local typ = H5.H5Dget_type(self._hid)
   local ret = hdf5.DataType(typ) -- a copy of typ is made
   H5.H5Tclose(typ)
   return ret
end

function hdf5.DataSet:set_mpio(mode)
   self._mpio.requested_mode = mode
end

function hdf5.DataSet:get_mpio()
   local ret = { }
   for k,v in pairs(self._mpio) do ret[k] = v end
   return ret
end

function hdf5.DataSet:__tostring__()
   if self._hid ~= 0 then
      return string.format("<HDF5 data set: \"%s\">", self._name)
   else
      return "<Closed HDF5 data set>"
   end
end

function hdf5.DataSet:__index__(slice)
   local slice = slice or { }
   if self._hid == 0 then
      return nil
   elseif type(slice) == 'table' then
      local htype = self:get_type()
      if htype:type_class() ~= 'float' then
	 error("DataSet:can only index data sets whose class is 'float'")
      end
      --------------------------------------------------------------------------
      -- A[slice] := A[{{i0,i1,di}, {j0:j1,dj}}] := A[i0:i1:di, j0:j1:dj]
      -- 
      -- defaults
      -- 
      -- i0 = i0 or 0
      -- i1 = i1 or Ni-1
      -- di = di or 1
      --------------------------------------------------------------------------
      local fspace = self:get_space()
      local mspace = hdf5.DataSpace()
      local extent = fspace:get_extent()
      local rank = #extent
      local start, stride, count, block = { }, { }, { }, { }
      for i=1,rank do
	 slice[i] = slice[i] or { }
	 start[i] = slice[i][1] or 0
	 stride[i] = slice[i][3] or 1
	 count[i] = ((slice[i][2] or extent[i]) - start[i]) / stride[i]
	 block[i] = 1
      end
      mspace:set_extent(count)
      fspace:select_hyperslab(start, stride, count, block)
      local buf = buffer.new_buffer(mspace:get_select_npoints() *
				    self:get_type():get_size())
      self:read(buf, mspace, fspace)
      return array.view(buf, htype:type_string(), count)
   end
end

function hdf5.DataSet:__newindex__(slice, value)
   --------------------------------------------------------------------------
   -- The variables `slice` and `value` must conform to:
   --
   -- slice = {{i0,i1,di}, {j0,j1,dj}}
   -- start, stride, count, block = value:selection()
   -- dtype = value:dtype()
   -- buffer = value:buffer()
   --
   --------------------------------------------------------------------------
   local slice = slice or { } -- whole selection
   if self._hid == 0 then
      error("DataSet:cannot write to closed data set")

   elseif type(slice) == 'string' then
      return Base.__index__(self, slice)

   elseif type(slice) == 'table' then
      local htype = self:get_type()
      local buf = value:buffer()

      if htype:type_class() ~= 'float' then
	 error("DataSet:can only assign to data sets with class 'float'")
      elseif htype:type_string() ~= value:dtype() then
	 error('DataSet:source and destination must have the same data type')
      end

      --------------------------------------------------------------------------
      -- Parse the `slice` table in order to determine target (file space)
      -- memory selection.
      --------------------------------------------------------------------------
      local fspace = self:get_space()
      local extent = fspace:get_extent()
      local rank = #extent
      local start, stride, count, block = { }, { }, { }, { }

      for i=1,rank do
	 slice[i] = slice[i] or { }
	 start[i] = slice[i][1] or 0
	 stride[i] = slice[i][3] or 1
	 count[i] = ((slice[i][2] or extent[i]) - start[i]) / stride[i]
	 block[i] = 1
      end
      fspace:select_hyperslab(start, stride, count, block)

      --------------------------------------------------------------------------
      -- Build the source memory selection from `value` object description.
      --------------------------------------------------------------------------
      local start, stride, count, block = value:selection()
      local mspace = hdf5.DataSpace(value:extent())
      mspace:select_hyperslab(start, stride, count, block)
      self:write(buf, mspace, fspace)
   else
      error('DataSet:index object not recognized')
   end
end

--------------------------------------------------------------------------------
-- HDF5 DataType class methods
--------------------------------------------------------------------------------
function hdf5.DataType:get_size()
   return H5.H5Tget_size(self._hid)
end

function hdf5.DataType:set_size(size)
   return H5.H5Tset_size(self._hid, size)
end

function hdf5.DataType:type_string()
   local htyp = self._hid
   local tcls = H5.H5Tget_class(htyp)
   if tcls == H5.H5T_STRING then return 'string'
   elseif H5.H5Tequal(htyp, H5.H5T_NATIVE_CHAR) then return 'char'
   elseif H5.H5Tequal(htyp, H5.H5T_NATIVE_INT) then return 'int'
   elseif H5.H5Tequal(htyp, H5.H5T_NATIVE_FLOAT) then return 'float'
   elseif H5.H5Tequal(htyp, H5.H5T_NATIVE_DOUBLE) then return 'double'
   else return nil
   end
end

function hdf5.DataType:type_class()
   return ({[H5.H5T_INTEGER]='int',
	    [H5.H5T_FLOAT]='float',
	    [H5.H5T_STRING]='string',
	    [H5.H5T_BITFIELD]='bitfield',
	    [H5.H5T_OPAQUE]='opaque',
	    [H5.H5T_COMPOUND]='compound',
	    [H5.H5T_REFERENCE]='reference',
	    [H5.H5T_ENUM]='enum',
	    [H5.H5T_VLEN]='vlen',
	    [H5.H5T_ARRAY]='array'})[H5.H5Tget_class(self._hid)]
end

function hdf5.DataType:__tostring__()
   if self._hid ~= 0 then
      return string.format("<HDF5 data type: \"%s\">", self._name)
   else
      return "<Closed HDF5 data type>"
   end
end


--------------------------------------------------------------------------------
-- HDF5 DataSpace class methods
--------------------------------------------------------------------------------
function hdf5.DataSpace:get_extent(what)
   if self._hid == 0 then
      error("DataSpace:get_extent cannot operate on closed data space")
   end
   local rank = H5.H5Sget_simple_extent_ndims(self._hid)
   local csize = { }
   local msize = { }
   for i=1,rank do
      csize[i] = 0
      msize[i] = 0
   end
   local current_size = H5.new_hsize_t_arr(csize)
   local maximum_size = H5.new_hsize_t_arr(msize)
   H5.H5Sget_simple_extent_dims(self._hid, current_size, maximum_size)
   for i=1,rank do
      csize[i] = current_size[i-1]
      msize[i] = maximum_size[i-1]
   end
   return ({current=csize, maximum=msize})[what] or csize
end

function hdf5.DataSpace:get_select_npoints()
   return H5.H5Sget_select_npoints(self._hid)
end

function hdf5.DataSpace:set_extent(extent, max)
   if self._hid == 0 then
      error("DataSpace:set_extent cannot operate on closed data space")
   end
   if max and #extent ~= #max then
      error("DataSpace:extent and max have different sizes")
   end
   local current_size = H5.new_hsize_t_arr(extent)
   local maximum_size = H5.new_hsize_t_arr(max or extent)
   local err = H5.H5Sset_extent_simple(self._hid, #extent, current_size,
				       maximum_size)
   if #err < 0 then error("DataSpace:set_extent") end
end

function hdf5.DataSpace:select_hyperslab(start, stride, count, block)
   if self._hid == 0 then
      error("DataSpace:select_hyperslab cannot operate on closed data space")
   end
   local hstart = H5.new_hsize_t_arr(start)
   local hstride = H5.new_hsize_t_arr(stride)
   local hcount = H5.new_hsize_t_arr(count)
   local hblock = H5.new_hsize_t_arr(block)
   local err = H5.H5Sselect_hyperslab(self._hid, H5.H5S_SELECT_SET, hstart, 
				      hstride, hcount, hblock)
   if #err < 0 then error("DataSpace:select_hyperslab") end
end

function hdf5.DataSpace:__tostring__()
   if self._hid ~= 0 then
      return string.format("<HDF5 data space: [%s]>",
			   table.concat(self:get_extent(), " "))
   else
      return "<Closed HDF5 data space>"
   end
end


--------------------------------------------------------------------------------
-- HDF5 File constructor
--------------------------------------------------------------------------------
function hdf5.File:__init__(name, mode, opts)
   oo.setattrib(self, '_parent', parent)
   oo.setattrib(self, '_name', name)
   oo.setattrib(self, '_type', 'file')
   oo.setattrib(self, '_hid', 0)
   oo.setattrib(self, '_close', H5.H5Fclose)
   oo.setattrib(self, '_open_objects', { })

   local opts = opts or { }
   local fcpl = H5.H5Pcreate(H5.H5P_FILE_CREATE)
   local fapl = H5.H5Pcreate(H5.H5P_FILE_ACCESS)

   if opts.mpi then
      H5.H5Pset_fapl_mpio(fapl, opts.mpi.comm, opts.mpi.info)
   end
   if opts.align then
      H5.H5Pset_alignment(fapl, opts.align.threshold, opts.align.alignment)
   end
   if opts.btree_ik then
      H5.H5Pset_istore_k(fcpl, opts.btree_ik)
   end

   if mode == "w" then
      self._hid = H5.H5Fcreate(name, H5.H5F_ACC_TRUNC, fcpl, fapl)
   elseif mode == "r" then
      self._hid = H5.H5Fopen(name, H5.H5F_ACC_RDONLY, fapl)
   elseif mode == "r+" then
      self._hid = H5.H5Fopen(name, H5.H5F_ACC_RDWR, fapl)
   else
      error("File:mode must be one of [w, r, r+]")
   end

   if #self._hid < 0 then error("File:could not open or create file") end
   H5.H5Pclose(fcpl)
   H5.H5Pclose(fapl)
end


--------------------------------------------------------------------------------
-- HDF5 Group constructor
--------------------------------------------------------------------------------
function hdf5.Group:__init__(parent, name)
   oo.setattrib(self, '_parent', parent)
   oo.setattrib(self, '_name', name)
   oo.setattrib(self, '_type', 'group')
   oo.setattrib(self, '_hid', 0)
   oo.setattrib(self, '_close', H5.H5Gclose)
   oo.setattrib(self, '_open_objects', { })

   if not H5.H5Lexists(parent._hid, name, hp0) then
      self._hid = H5.H5Gcreate2(parent._hid, name, hp0, hp0, hp0)
   else
      self._hid = H5.H5Gopen2(parent._hid, name, hp0)
   end
   parent._open_objects[name] = self
end


--------------------------------------------------------------------------------
-- HDF5 DataSet constructor
--------------------------------------------------------------------------------
function hdf5.DataSet:__init__(parent, name, mode, opts)
   oo.setattrib(self, '_parent', parent)
   oo.setattrib(self, '_name', name)
   oo.setattrib(self, '_type', 'data set')
   oo.setattrib(self, '_hid', 0)
   oo.setattrib(self, '_mpio', { })
   oo.setattrib(self, '_close', H5.H5Dclose)

   local opts = opts or { }
   local mode = mode or "r+"

   if mode == "w" then
      local space = hdf5.DataSpace(opts.shape, opts.max)
      local dtype = hdf5.DataType(opts.dtype)

      if H5.H5Lexists(parent._hid, name, hp0) then
       	 local err = H5.H5Ldelete(parent._hid, name, hp0)
       	 if #err < 0 then
       	    error("DataSet:failed to clobber existing data set")
       	 else
       	    print("DataSet:successfully clobbered existing data set")
       	 end
      end
      local dcpl = H5.H5Pcreate(H5.H5P_DATASET_CREATE)
      if opts.chunk then
       	 local c = H5.new_hsize_t_arr(opts.chunk)
       	 local err = H5.H5Pset_chunk(dcpl, #opts.chunk, c)
	 if #err < 0 then error("DataSet:could not set chunk") end
      end
      self._hid = H5.H5Dcreate2(parent._hid, name, dtype._hid, space._hid, hp0,
      			       dcpl, hp0)
      H5.H5Pclose(dcpl)
   elseif mode == "r+" then
      if not H5.H5Lexists(parent._hid, name, hp0) then
	 error("DataSet:cannot open data set for reading")
      end
      self._hid = H5.H5Dopen2(parent._hid, name, hp0)
   else
      error("DataSet:mode must be one of [w, r+]")
   end

   if #self._hid < 0 then
      error("DataSet:error opening or creating data set")
   end
   parent._open_objects[name] = self
end


--------------------------------------------------------------------------------
-- HDF5 DataType constructor
--
-- If `arg` is an hid_t representing a type then a copy of that type is made via
-- H5Tcopy. If `arg` is a string, then it must be a valid key into the variable
-- typedict. Opening existing types from data sets is done by the DataSet class.
--------------------------------------------------------------------------------
function hdf5.DataType:__init__(typeid)
   oo.setattrib(self, '_type', 'data type')
   oo.setattrib(self, '_hid', 0)
   oo.setattrib(self, '_close', H5.H5Tclose)

   local typedict = {char=H5.H5T_NATIVE_CHAR,
		     int=H5.H5T_NATIVE_INT,
		     float=H5.H5T_NATIVE_FLOAT,
		     double=H5.H5T_NATIVE_DOUBLE}
   if type(typeid) == 'string' then
      self._hid = H5.H5Tcopy(typedict[typeid])
   else
      self._hid = H5.H5Tcopy(typeid)
   end
end


--------------------------------------------------------------------------------
-- HDF5 DataSpace constructor
--
-- If `arg` is a string, then it must be either 'simple' or 'scalar' and a new
-- data space is returned [default=simple]. If `arg` is a numeric table then a
-- new simple data space with that extent (and maximum extent `max` if present)
-- is returned. If `arg` is a data set then its data space is returned.
--------------------------------------------------------------------------------
function hdf5.DataSpace:__init__(arg, max)
   oo.setattrib(self, '_type', 'data space')
   oo.setattrib(self, '_hid', 0)
   oo.setattrib(self, '_close', H5.H5Sclose)

   if not arg or type(arg) == 'string' then
      local t = { simple=H5.H5S_SIMPLE, scalar=H5.H5S_SCALAR }
      self._hid = H5.H5Screate(t[arg or 'simple'])
   elseif arg._type == 'data set' then
      self._hid = H5.H5Dget_space(arg._hid)
   elseif type(arg) == 'table' then
      self._hid = H5.H5Screate(H5.H5S_SIMPLE)
      self:set_extent(arg, max)
   else
      error("DataSpace:constructor argument not understood")
   end
   if #self._hid < 0 then error("DataSpace:creation failed") end
end


--------------------------------------------------------------------------------
-- Unit test cases
--------------------------------------------------------------------------------
local function test1()
   local h5f = hdf5.File("outfile.h5", "w")
   local h5g = hdf5.Group(h5f, "thegroup")
   local h5h = hdf5.Group(h5g, "thesubgroup")
   assert(h5f["thegroup"]["thesubgroup"] == h5h)
   h5f:close()
end

local function test2()
   local h5f = hdf5.File("outfile.h5", "w")
   local h5g = hdf5.Group(h5f, "thegroup")
   local h5h = hdf5.Group(h5g, "thesubgroup")
   h5g:close()
   assert(h5f["thegroup"]["thesubgroup"]:path() == h5h:path())
   assert(h5f["thething"] == nil)
   h5f:close()
end

local function test3()
   local h5f = hdf5.File("outfile.h5", "w")
   local hg1 = hdf5.Group(h5f, "thegroup1")
   local hg2 = hdf5.Group(h5f, "thegroup2")
   local hg3 = hdf5.Group(h5f, "thegroup3")
   for k,v in pairs(h5f) do
      assert(v == h5f[k])
   end
   hg3["message"] = "the message content"
   h5f:close()
end

local function test4()
   local h5s = hdf5.DataSpace()
   h5s:set_extent{10,12,14}
   local size = h5s:get_extent('maximum')
   assert(size[1] == 10)
   assert(size[2] == 12)
   assert(size[3] == 14)

   local dtype = hdf5.DataType('double')
   dtype:set_size(128)
   assert(dtype:get_size() == 128)
end

local function test5() -- depends on test3 being run first
   local file = hdf5.File("outfile.h5", "r")
   local dset = hdf5.DataSet(file["thegroup3"], "message")
   local fspc = hdf5.DataSpace(dset)
   local buf = dset:read()
   assert(tostring(buf) == "the message content")
   assert(type(file["thegroup3"]["message"]:value()) == 'string')
   dset:close()
   file:close()
end

local function test6()
   local file = hdf5.File("outfile.h5", "w")
   local dset = hdf5.DataSet(file, "data1d", "w",
			     {shape={4,4,8}, dtype='double'})
   local buf = buffer.new_buffer(4*4*8*8)
   dset:write(buf)

   local array = array.view(buf, 'double', {4,4,8}, {0,0,0}, {2,2,2}, {2,2,2})
   file["data3d"] = array
   file:close()

   local file = hdf5.File("outfile.h5", "r")
   assert(#file["data3d"]:read() == 64)
   assert(file["data3d"]:value():dtype() == 'double')
   local start, stride, count, block = file["data3d"]:value():selection()
   assert(count[1] == 2)
   assert(count[2] == 2)
   assert(count[3] == 2)
   file:close()
end

local function test7()
   local buf = buffer.new_buffer(4*4*8 * array.sizeof('double'))
   local my_data = array.view(buf, 'double', {4,4,8})
   local h5f = hdf5.File("outfile.h5", "w")
   h5f["dataset"] = my_data
   local group1 = h5f:require_group("group1")
   group1["message"] = "here is the message"

   assert(group1["message"]:get_type():type_class() == 'string')
   assert(group1["message"]:get_type():type_string() == 'string')
   assert(h5f["dataset"]:get_type():type_class() == 'float')
   assert(h5f["dataset"]:get_type():type_string() == 'double')

   local space = hdf5.DataSpace()
   space:set_extent{4,4,8}
   space:select_hyperslab({0,0,0}, {1,1,1}, {4,4,8}, {1,1,1})
   h5f["dataset"]:read(buf, space, space)
   local read_select = h5f["dataset"][{{0,4,2},{0,4,2},{0,8,2}}]
   assert(#read_select == 16)

   h5f["dataset"][{{0,4,2},{0,4,2},{0,4,1}}] = read_select
   local read_shape = read_select:shape()
   assert(read_shape[1] == 2)
   assert(read_shape[2] == 2)
   assert(read_shape[3] == 4)
   h5f:close()
end

local function test8()
   local h5f = hdf5.File("outfile.h5", "w")
   local h5d = hdf5.DataSet(h5f, "dataset", 'w',
			    {dtype='double', shape={10,10}, max={20,H5.H5S_UNLIMITED},
			     chunk={5,10}})
   h5f["dataset"]:set_extent{20,20}
   h5f:close()
   local h5f = hdf5.File("outfile.h5", "r")
   local chunk = h5f["dataset"]:get_chunk()
   assert(chunk[1] == 5)
   assert(chunk[2] == 10)
   h5f:close()

   local h5f = hdf5.File("outfile.h5", "w")
   local h5d = hdf5.DataSet(h5f, "dataset", 'w', {dtype='double', shape={10,10}})
   h5f:close()
   local h5f = hdf5.File("outfile.h5", "r")
   local chunk = h5f["dataset"]:get_chunk()
   assert(#chunk == 0)
   h5f:close()
end


if ... then -- if __name__ == "__main__"
   return hdf5
else
   test1()
   test2()
   test3()
   test4()
   test5()
   test6()
   test7()
   test8()
   print(debug.getinfo(1).source, ": All tests passed")
end
