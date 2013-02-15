
local hdf5   = require 'lua-hdf5.LuaHDF5'
local array  = require 'array'
local visual = require 'visual'
local util   = require 'util'

for _,fname in pairs(arg) do
   if util.endswith(fname, '.h5') then
      local h5f = hdf5.File(fname, 'r')
      local D = h5f['prim']['rho']:value()
      local N = D:shape()
      local ppmname = string.gsub(fname, '.h5', '.ppm')
      local pngname = string.gsub(fname, '.h5', '.png')
      print('writing '..ppmname)
      visual.write_ppm(D:buffer(), N[1], N[2], ppmname)
      os.execute(string.format('convert %s %s', ppmname, pngname))
      os.execute(string.format('rm %s', ppmname))
   end
end
