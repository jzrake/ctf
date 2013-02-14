
local oo   = require 'class'
local json = require 'json'
local util = { }


-- *****************************************************************************
-- Function to deep-copy a table
-- .............................................................................
function util.deepcopy(object)
   local lookup_table = {}
   local function _copy(object)
      if type(object) ~= "table" then
         return object
      elseif lookup_table[object] then
         return lookup_table[object]
      end
      local new_table = {}
      lookup_table[object] = new_table
      for index, value in pairs(object) do
         new_table[_copy(index)] = _copy(value)
      end
      return setmetatable(new_table, getmetatable(object))
   end
   return _copy(object)
end


-- *****************************************************************************
-- Check if a string `s` end swith the string `send`
-- .............................................................................
function util.endswith(s, send)
   return #s >= #send and s:find(send, #s-#send+1, true) and true or false
end


-- *****************************************************************************
-- Utility to print tables
-- .............................................................................
function util.pretty_print(t, indent, compare)
   local names = { }
   if not indent then indent = "" end
   for n,g in pairs(t) do
      table.insert(names,n)
   end
   table.sort(names, compare)
   for i,n in pairs(names) do
      local v = t[n]
      if oo.isclass(v) then
	 print(indent..oo.classname(v))
      elseif type(v) == "table" then
	 if(v==t) then
	    print(indent..tostring(n)..": <-")
	 else
	    print(indent..tostring(n)..":")
	    util.pretty_print(v,indent.."   ")
	 end
      else
	 if type(v) == "function" then
	    print(indent..tostring(n).."()")
	 elseif type(v) == 'number' and (math.log10(v) > 4 or
					 math.log10(v) < -4) then
	    print(indent..tostring(n)..": "..string.format("%3.2e", v))
	 else
	    print(indent..tostring(n)..": "..tostring(v))
	 end
      end
   end
end


-- *****************************************************************************
-- Function to call Gnuplot from Lua using popen
-- .............................................................................
function util.plot(series, opts)
   local opts = opts or { }
   local gp = io.popen("gnuplot", 'w')

   if opts.output then
      gp:write("set terminal postscript enhanced color\n")
      gp:write(string.format("set output '| ps2pdf - %s.pdf'\n", opts.output))
      gp:write(string.format("set title '%s'\n", opts.id))
   end

   for _,cmd in pairs(opts.cmds or { }) do
      gp:write(cmd..'\n')
   end
   local lines = { }
   for k,v in pairs(series) do
      table.insert(lines, string.format(" '-' u 1:2 %s title '%s'",
					opts.ls or 'w lp', k))
   end
   gp:write("plot" .. table.concat(lines, ",") .. "\n")
   for k,v in pairs(series) do
      if #v == 0 then
	 local xy_vec = { }
	 for x,y in pairs(v) do
	    table.insert(xy_vec, {x,y})
	 end
	 table.sort(xy_vec, function(a,b) return a[1] < b[1] end)
	 for _,pnt in pairs(xy_vec) do
	    gp:write(string.format("%12.10e %12.10e\n", pnt[1], pnt[2]))
	 end
      else
	 for x,y in pairs(v) do
	    gp:write(string.format("%12.10e %12.10e\n", x, y))
	 end
      end
      gp:write("e\n")
   end
   if not opts.output then
      gp:write(string.format(" pause %f\n", opts.tpause or 100.0))
   end
   gp:close()
end


-- *****************************************************************************
-- Function to parse cmdline.opts into the runargs table
-- .............................................................................
function util.parse_args(runargs)
   for k,v in pairs(cmdline.opts) do
      if type(runargs[k]) == 'number' then
         runargs[k] = tonumber(v)
      elseif type(runargs[k]) == 'boolean' then
	 runargs[k] = tonumber(v) == 1 or v == "true"
      else
         runargs[k] = v
      end
   end
   print("command line arguments:")
   for k,v in pairs(runargs) do
      print(k, " ........ ", v)
   end
end


-- *****************************************************************************
-- http://lua-users.org/wiki/SplitJoin
-- -----------------------------------------------------------------------------
function util.string_split(self, sSeparator, nMax, bRegexp)
   assert(sSeparator ~= '')
   assert(nMax == nil or nMax >= 1)

   local aRecord = {}

   if self:len() > 0 then
      local bPlain = not bRegexp
      nMax = nMax or -1

      local nField=1 nStart=1
      local nFirst,nLast = self:find(sSeparator, nStart, bPlain)
      while nFirst and nMax ~= 0 do
         aRecord[nField] = self:sub(nStart, nFirst-1)
         nField = nField+1
         nStart = nLast+1
         nFirst,nLast = self:find(sSeparator, nStart, bPlain)
         nMax = nMax-1
      end
      aRecord[nField] = self:sub(nStart)
   end

   return aRecord
end

return util
