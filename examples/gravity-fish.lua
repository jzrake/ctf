
local array  = require 'array'
local fish   = require 'fish'
local fluids = require 'fluids'
local hdf5   = require 'lua-hdf5.LuaHDF5'

local P = array.array{100, 5}
local G = array.array{100, 4}
fish.run_euler(P:buffer(), G:buffer())

local outfile = hdf5.File('euler.h5', 'w')
outfile['prim'] = P
outfile['grav'] = G
outfile:close()

os.execute("python plot.py")
