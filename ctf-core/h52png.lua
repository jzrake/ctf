
if arg[2] == '--explain' then
   print "Convert 2d HDF5 data to images, optionally make a movie with ffmpeg"
   return
end

local hdf5     = require 'lua-hdf5.LuaHDF5'
local array    = require 'array'
local visual   = require 'visual'
local util     = require 'util'
local optparse = require 'optparse'

local frames_written = { }

local function file_exists(name)
   local f = io.open(name,"r")
   if f ~= nil then io.close(f) return true else return false end
end

local function write_frames(opts, fnames)
   for _,fname in pairs(fnames) do
      local h5f = hdf5.File(fname, 'r')
      local D = h5f['prim']['rho']:value()
      local N = D:shape()
      local ppmname = string.gsub(fname, '.h5', '.ppm')
      local pngname = string.gsub(fname, '.h5', '.png')
      print('[h52png]: '..fname..' -> '..pngname)
      visual.write_ppm(D:buffer(), N[1], N[2], ppmname,
                       opts.cmap, opts.dmin, opts.dmax)
      os.execute(string.format('convert -flip %s %s', ppmname, pngname))
      os.execute(string.format('rm %s', ppmname))
      table.insert(frames_written, pngname)
   end
end

local function main()
   local usage = "h52png <problem> [<options>]"
   local parser = optparse.OptionParser{usage=usage,
                                        version="CTF version 1.0"}
   parser.add_option{"--dmin", dest="dmin", help="data range minimum"}
   parser.add_option{"--dmax", dest="dmax", help="data range maximum"}
   parser.add_option{"--cmap", dest="cmap", help="color map index"}
   parser.add_option{"--movie", dest="movie", help="write a movie with the given name"}
   parser.add_option{"--format", dest="format",
                     help="use instead of direct filenames: for example, data/myrun.%04d"}
   parser.add_option{"--cleanup", dest="cleanup", help="remove the images after finished"}

   local opts, args = parser.parse_args()
   local fnames = { }
   if opts.format then
      local i=0
      while true do
         local testf = string.format(opts.format..'.h5', i)
         if file_exists(testf) then
            table.insert(fnames, testf)
            i = i + 1
         else
            break
         end
      end
   else
      for _,f in ipairs(args) do
         if util.endswith(f, '.h5') then
            table.insert(fnames, f)
         end
      end
   end
   write_frames(opts, fnames)

   if opts.movie then
      if not opts.format then
	 error("format option required: something like --format=data/myrun.%04d")
      end
      if not util.endswith(opts.movie, '.mp4') then
	 error("the movie extension should be mp4")
      end
      os.execute(string.format("ffmpeg -f image2 -i %s.png -vcodec mpeg4 -b:v 1600k %s",
			       opts.format, opts.movie))
   end

   if opts.cleanup then
      for _,f in ipairs(frames_written) do
         os.execute(string.format('rm %s', f))
      end
   end
end

main()
