
local MPI    = require 'MPI'
local array  = require 'array'
local oo     = require 'class'
local cow    = require 'cow'
local hdf5   = require 'lua-hdf5.LuaHDF5'
local util   = require 'util'

local DataManagerHDF5 = oo.class('DataManagerHDF5')
function DataManagerHDF5:__init__(domain, dataset_names, opts)

   local opts = opts or { }
   local Nx = cow.domain_getsize(domain, 0)
   local Ny = cow.domain_getsize(domain, 1)
   local Nz = cow.domain_getsize(domain, 2)
   local Ng = cow.domain_getguard(domain)
   local i0 = cow.domain_getglobalstartindex(domain, 0)
   local j0 = cow.domain_getglobalstartindex(domain, 1)
   local k0 = cow.domain_getglobalstartindex(domain, 2)
   local i1 = cow.domain_getnumlocalzonesinterior(domain, 0) + i0
   local j1 = cow.domain_getnumlocalzonesinterior(domain, 1) + j0
   local k1 = cow.domain_getnumlocalzonesinterior(domain, 2) + k0
   local Ni = i1 - i0
   local Nj = j1 - j0
   local Nk = k1 - k0
   local Nq = #dataset_names

   self.ndim = cow.domain_getndim(domain)
   self.domain = domain
   self.dataset_names = dataset_names
   self.Ng = Ng

   if self.ndim == 1 then
      self.array_shape = {Ni+2*Ng, #dataset_names}
      self.ggrid_shape = {Nx} -- global grid shape
      self.sgrid_shape = {Ni} -- subgrid shape
      self.sgrid_start = {i0}
      self.array = array.array(self.array_shape, 'double')
   elseif self.ndim == 2 then
      self.array_shape = {Ni+2*Ng, Nj+2*Ng, #dataset_names}
      self.ggrid_shape = {Nx, Ny} -- global grid shape
      self.sgrid_shape = {Ni, Nj} -- subgrid shape
      self.sgrid_start = {i0, j0}
      self.array = array.array(self.array_shape, 'double')
   elseif self.ndim == 3 then
      self.array_shape = {Ni+2*Ng, Nj+2*Ng, Nk+2*Ng, #dataset_names}
      self.ggrid_shape = {Nx, Ny, Nz} -- global grid shape
      self.sgrid_shape = {Ni, Nj, Nk} -- subgrid shape
      self.sgrid_start = {i0, j0, k0}
      self.array = array.array(self.array_shape, 'double')
   end
   local domain_comm = MPI.Comm()
   cow.domain_getcomm(domain, domain_comm)

   -----------------------------------------------------------------------------
   -- Configuration options: try many to fine-tune performance
   -----------------------------------------------------------------------------
   self.file_opts = { }
   self.file_opts.mpi = {comm = domain_comm,
			 info = MPI.INFO_NULL}
   self.file_opts.align = {threshold = opts.align and 4 * KB or 1,
			   alignment = opts.align and stripe_size or 1}
   self.file_opts.btree_ik = opts.btree_ik or 32 -- default

   self.dset_opts = { }
   self.dset_opts.chunk = opts.chunk or false
   self.dset_opts.mpio = opts.mpio
end

function DataManagerHDF5:synchronize()
   local field = cow.dfield_new()
   for _,v in ipairs(self.dataset_names) do
      cow.dfield_addmember(field, v)
   end
   cow.dfield_setdomain(field, self.domain)
   cow.dfield_setdatabuffer(field, self.array:buffer())
   cow.dfield_commit(field)
   cow.dfield_syncguard(field)
   cow.dfield_del(field)
end

function DataManagerHDF5:local_mesh_size(which)
   local which = which or 'total'
   if which == 'total' then
      return cow.domain_getnumlocalzonesinterior(self.domain, cow.ALL_DIMS)
   elseif which == 'shape' then
      local size = { }   
      for i=1,self.ndim do
	 size[i] = cow.domain_getnumlocalzonesinterior(self.domain, i-1)
      end
      return size
   else
      error("argument 'which' must be ['total', 'shape']")
   end
end

function DataManagerHDF5:local_extent(which)
   local lower = { }
   local upper = { }
   for i=1,self.ndim do
      lower[i] = cow.domain_getlowercoord(self.domain, i-1)
      upper[i] = cow.domain_getuppercoord(self.domain, i-1)
   end
   return lower, upper
end

function DataManagerHDF5:write(filename, opts)
   if not self.dset_opts.mpio then
      return self:write_sequential(filename, opts)
   else
      return self:write_parallel(filename, opts)
   end
end

function DataManagerHDF5:write_parallel(filename, opts)
   local opts = opts or { }
   local start = os.clock()
   local file = hdf5.File(filename, opts.file_mode or 'w', self.file_opts)
   local group = opts.group and hdf5.Group(file, opts.group) or file
   local dset_opts = {shape=self.ggrid_shape, dtype='double'}

   if self.dset_opts.chunk then
      dset_opts.chunk = self.sgrid_shape
   end

   print("[HDF5-parallel] writing to file " .. filename)

   for i,name in ipairs(self.dataset_names) do
      print("[HDF5-parallel] writing data set " .. name)
      local dset = hdf5.DataSet(group, name, opts.dset_mode or 'w', dset_opts)
      local mspace = hdf5.DataSpace(self.array_shape)
      local fspace = dset:get_space()

      self:_setup_spaces(mspace, fspace, i)
      dset:set_mpio(self.dset_opts.mpio)
      dset:write(self.array:buffer(), mspace, fspace)

      if i == #self.dataset_names  then
	 print("[HDF5-parallel] write stats:")
	 util.pretty_print(dset:get_mpio(), '\t# ')
      end
      dset:close()
   end

   file:close()
   MPI.Barrier(self.file_opts.mpi.comm)

   local dt = os.clock() - start
   print(string.format("[HDF5-parallel] write time: %3.2f seconds", dt))

   return dt
end

function DataManagerHDF5:write_sequential(filename, opts)
   local start = os.clock()
   local size = cow.domain_getcartsize(self.domain)
   local rank = cow.domain_getcartrank(self.domain)
   local dset_opts = {shape=self.ggrid_shape, dtype='double'}

   if self.dset_opts.chunk then
      dset_opts.chunk = self.sgrid_shape
   end

   if rank == 0 and (opts.file_mode or 'w') == 'w' then
      local file = hdf5.File(filename, 'w')
      if (opts.dset_mode or 'w') == 'w' then
	 local group = opts.group and hdf5.Group(file, opts.group) or file
	 for i,name in ipairs(self.dataset_names) do
	    local dset = hdf5.DataSet(group, name, 'w', dset_opts)
	 end
      end
      file:close()
   end

   for irank=0,size-1 do
      if irank == rank then
   	 local file = hdf5.File(filename, 'r+')
   	 local group = opts.group and hdf5.Group(file, opts.group) or file

   	 print("[HDF5-sequential] writing to file " .. filename)
	 
   	 for i,name in ipairs(self.dataset_names) do
   	    print("[HDF5-sequential] writing data set " .. name)
   	    local dset = hdf5.DataSet(group, name, 'r+', dset_opts)
   	    local mspace = hdf5.DataSpace(self.array_shape)
   	    local fspace = dset:get_space()
	    
   	    self:_setup_spaces(mspace, fspace, i)
   	    dset:write(self.array:buffer(), mspace, fspace)
   	    dset:close()
   	 end
   	 file:close()
      end
      MPI.Barrier(self.file_opts.mpi.comm)
   end

   local dt = os.clock() - start
   print(string.format("[HDF5-sequential] write time: %3.2f seconds", dt))

   return dt
end

function DataManagerHDF5:read(filename, opts)
   local opts = opts or { }
   local start = os.clock()
   local file = hdf5.File(filename, 'r', self.file_opts)
   local group = opts.group and hdf5.Group(file, opts.group) or file

   print("[HDF5] reading from file " .. filename)

   for i,name in ipairs(self.dataset_names) do
      print("[HDF5] reading data set " .. name)
      local dset = hdf5.DataSet(group, name, 'r+')
      local mspace = hdf5.DataSpace(self.array_shape)
      local fspace = dset:get_space()

      self:_setup_spaces(mspace, fspace, i)
      dset:set_mpio(self.dset_opts.mpio)
      dset:read(self.array:buffer(), mspace, fspace)

      if i == #self.dataset_names then
	 print("[HDF5] read stats:")
	 util.pretty_print(dset:get_mpio(), '\t# ')
      end
      dset:close()
   end

   file:close()
   MPI.Barrier(self.file_opts.mpi.comm)

   local dt = os.clock() - start
   print(string.format("[HDF5] read time: %3.2f seconds", dt))
   return dt
end

function DataManagerHDF5:power_spectrum(nbins)
   if #self.dataset_names ~= 3 then
      error("[HDF5] need a three-dimensional field to get a power spectrum")
   end

   local nbins = nbins or 128
   local binloc = array.vector(nbins)
   local binval = array.vector(nbins)
   local field = cow.dfield_new()
   local pspec = cow.histogram_new()

   for _,v in ipairs(self.dataset_names) do
      cow.dfield_addmember(field, v)
   end

   cow.dfield_setdomain(field, self.domain)
   cow.dfield_setdatabuffer(field, self.array:buffer())
   cow.dfield_commit(field)
   cow.histogram_setnbins(pspec, 0, nbins)
   cow.histogram_setspacing(pspec, cow.HIST_SPACING_LINEAR)
   cow.fft_pspecvecfield(field, pspec)
   cow.histogram_getbinlocx(pspec, binloc:buffer())
   cow.histogram_getbinvalv(pspec, binval:buffer())
   cow.histogram_del(pspec)
   cow.dfield_del(field)

   return binloc, binval
end

function DataManagerHDF5:_setup_spaces(mspace, fspace, member)
   local Ng = self.Ng
   local S = self.sgrid_shape
   local T = self.sgrid_start
   local Ni, Nj, Nk = S[1], S[2], S[3]
   local i0, j0, k0 = T[1], T[2], T[3]

   if self.ndim == 1 then
      local mstart = {Ng, member-1}
      local mstrid = {1, #self.dataset_names}
      local mcount = {Ni, 1}
      local mblock = {1, 1}
      local fstart = {i0}
      local fstrid = {1}
      local fcount = {Ni}
      local fblock = {1}
      mspace:select_hyperslab(mstart, mstrid, mcount, mblock)
      fspace:select_hyperslab(fstart, fstrid, fcount, fblock)
   elseif self.ndim == 2 then
      local mstart = {Ng, Ng, member-1}
      local mstrid = {1, 1, #self.dataset_names}
      local mcount = {Ni, Nj, 1}
      local mblock = {1, 1, 1}
      local fstart = {i0, j0}
      local fstrid = {1, 1}
      local fcount = {Ni, Nj}
      local fblock = {1, 1}
      mspace:select_hyperslab(mstart, mstrid, mcount, mblock)
      fspace:select_hyperslab(fstart, fstrid, fcount, fblock)
   elseif self.ndim == 3 then
      local mstart = {Ng, Ng, Ng, member-1}
      local mstrid = {1, 1, 1, #self.dataset_names}
      local mcount = {Ni, Nj, Nk, 1}
      local mblock = {1, 1, 1, 1}
      local fstart = {i0, j0, k0}
      local fstrid = {1, 1, 1}
      local fcount = {Ni, Nj, Nk}
      local fblock = {1, 1, 1}
      mspace:select_hyperslab(mstart, mstrid, mcount, mblock)
      fspace:select_hyperslab(fstart, fstrid, fcount, fblock)
   end
end
return {DataManagerHDF5=DataManagerHDF5}
