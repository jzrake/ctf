

-- *****************************************************************************
-- Function to deep-copy a table
-- .............................................................................
local function deepcopy(object)
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
-- Function to call Gnuplot from Lua using popen
-- .............................................................................
local function plot(series, tpause)
   local gp = io.popen("gnuplot", 'w')

   local lines = { }
   for k,v in pairs(series) do
      table.insert(lines, string.format(" '-' u 1:2 title '%s'", k))
   end
   
   gp:write("plot" .. table.concat(lines, ",") .. "\n")
   for k,v in pairs(series) do
      for i=0,v:size()-1 do
	 gp:write(string.format("%f %f\n", i, v[i]))
      end
      gp:write("e\n")
   end

   gp:write(string.format("pause %f\n", tpause or 100.0))
   gp:close()
end


return { deepcopy=deepcopy, plot=plot }
