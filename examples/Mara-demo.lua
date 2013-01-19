
local array   = require 'array'
local MPI     = require 'MPI'
local cow     = require 'cow'
local Mara    = require 'Mara'
local LuaMara = require 'Mara.LuaMara'
local hdf5    = require 'lua-hdf5.LuaHDF5'

MPI.Init()
cow.init(0, nil, 0) -- to reopen stdout to dev/null

Mara.start()
Mara.set_fluid('euler')
Mara.set_advance('single')
Mara.set_godunov('plm-muscl')
Mara.set_boundary('periodic')
Mara.set_riemann('hllc')

-- Global variables
local prim_names = Mara.fluid.GetPrimNames()
local Nx = 32
local Ny = 32
local Nz = 32
local Ng = 3
local Nq = #prim_names

local domain = cow.domain_new()
local domain_comm = MPI.Comm()
cow.domain_setndim(domain, 3)
cow.domain_setsize(domain, 0, Nx)
cow.domain_setsize(domain, 1, Ny)
cow.domain_setsize(domain, 2, Nz)
cow.domain_setguard(domain, Ng)
cow.domain_commit(domain)
cow.domain_getcomm(domain, domain_comm)

local function pinit(x,y,z)
   return {1, 1, math.sin(x), math.cos(y), 0}
end

local primitive = LuaMara.MaraDataManager(domain, prim_names)
local velocity = LuaMara.MaraDataManager(domain, {'vx', 'vy', 'vz'})
local P = primitive.array
local V = velocity.array

Mara.set_domain({0,0,0}, {1,1,1}, {Nx, Ny, Nz}, Nq, Ng, domain_comm)
Mara.init_prim(P:buffer(), pinit)
Mara.units.Print()

V[{nil,nil,nil,{0,1}}] = P[{nil,nil,nil,{2,3}}]
V[{nil,nil,nil,{1,2}}] = P[{nil,nil,nil,{3,4}}]
V[{nil,nil,nil,{2,3}}] = P[{nil,nil,nil,{4,5}}]

primitive:write('chkpt.0001.h5', {group='prim'})
velocity:write('velocity.h5')

local binloc, binval = velocity:power_spectrum(64)

if cow.domain_getcartrank(domain) == 0 then
   local h5f = hdf5.File('power_spectrum.h5', 'w')
   local grp = hdf5.Group(h5f, '0000')
   grp['binloc'] = binloc
   grp['binval'] = binval
   h5f:close()
end

--local time, error = Mara.advance(P:buffer(), 0.1)

cow.domain_del(domain)
Mara.close()
MPI.Finalize()
