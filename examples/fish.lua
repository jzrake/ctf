
local array  = require 'array'
local fish   = require 'fish'
local fluids = require 'fluids'
local hdf5   = require 'lua-hdf5.LuaHDF5'

local P = array.array{100, 5}
fish.run_euler(P:buffer())

local outfile = hdf5.File('euler.h5', 'w')
outfile['prim'] = P
outfile:close()

os.execute("python plot.py")
