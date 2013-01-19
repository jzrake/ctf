
--------------------------------------------------------------------------------
-- Test for the quality of the Poisson solver. The output is the relative L2
-- error, which should be less than or around 1%. The input is an HDF5 file
-- followed by the name of a three-dimensional data set from which to read a
-- source term field, rho. The equation
--
--                               del^2 phi = rho
--
-- is then solved solved for the potenial, phi using the FFT's in the cow
-- module. It is assumed (in the cow.transform.laplacian function) that the grid
-- spacing is equal along each axis, so make sure the domain physical size has
-- the same aspect ratio as the input data array. Periodic boundaries are
-- assumed on all walls.
--------------------------------------------------------------------------------
local array = require 'array'
local cow   = require 'cow'
local hdf5  = require 'lua-hdf5.LuaHDF5'

local fname = arg[2]
local dset = arg[3]

if not (fname or dset) then
   print "usage: poisson infile.h5 dataset_name"
   os.exit()
end

local infile = hdf5.File(fname, 'r')
local domain_shape = infile[dset]:get_space():get_extent()

local ng = 1 -- one ghost zone for laplacian
local domain = cow.domain_new()
cow.domain_setndim(domain, 3)
cow.domain_setsize(domain, 0, domain_shape[1])
cow.domain_setsize(domain, 1, domain_shape[2])
cow.domain_setsize(domain, 2, domain_shape[3])
cow.domain_setguard(domain, ng)
cow.domain_commit(domain)

local array_shape = { cow.domain_getnumlocalzonesincguard(domain, 0),
		      cow.domain_getnumlocalzonesincguard(domain, 1),
		      cow.domain_getnumlocalzonesincguard(domain, 2) }

local rho = cow.dfield_new()
local phi = cow.dfield_new()
local lph = cow.dfield_new()
local rho_data = array.array(array_shape)
local phi_data = array.array(array_shape)
local lph_data = array.array(array_shape) -- laplacian phi

rho_data[{{ng,-ng},{ng,-ng},{ng,-ng}}] = infile[dset][nil]
infile:close()

cow.dfield_setdomain(rho, domain)
cow.dfield_setdomain(phi, domain)
cow.dfield_setdomain(lph, domain)
cow.dfield_addmember(rho, "rho")
cow.dfield_addmember(phi, "phi")
cow.dfield_addmember(lph, "lph")
cow.dfield_setdatabuffer(rho, rho_data:buffer())
cow.dfield_setdatabuffer(phi, phi_data:buffer())
cow.dfield_setdatabuffer(lph, lph_data:buffer())
cow.dfield_commit(rho)
cow.dfield_commit(phi)
cow.dfield_commit(lph)

cow.dfield_syncguard(rho)
cow.fft_solvepoisson(rho, phi)
cow.dfield_syncguard(phi)

cow.dfield_clearargs(lph)
cow.dfield_pusharg(lph, phi)
cow.dfield_settransform(lph, cow.transform.laplacian)
cow.dfield_setuserdata(lph, nil)
cow.dfield_transformexecute(lph)
cow.dfield_syncguard(lph)

local rho_reduce = array.vector(3, 'double')
cow.dfield_clearargs(rho)
cow.dfield_settransform(rho, cow.transform.component)
cow.dfield_setuserdata(rho, cow.dfield_light(rho))
cow.dfield_setiparam(rho, 0)
cow.dfield_reduce(rho, rho_reduce:buffer())

local ninterior = cow.domain_getnumglobalzones(domain, cow.ALL_DIMS)
local rho_mean = rho_reduce[2] / ninterior

local rho_vec = rho_data:vector()
local lph_vec = lph_data:vector()
local dx = cow.domain_getgridspacing(domain, 0)
local L2 = 0.0
local VR = 0.0 -- variance in rho

for i=0,#rho_vec-1 do
   L2 = L2 + ((rho_vec[i] - rho_mean) - lph_vec[i]/dx^2)^2
   VR = VR + (rho_vec[i] - rho_mean)^2
end

print("relative L2 error:", L2 / VR)

cow.dfield_del(rho)
cow.dfield_del(phi)
cow.dfield_del(lph)
cow.domain_del(domain)
