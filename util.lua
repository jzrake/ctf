

local host = require 'host'

-- *****************************************************************************
-- Main driver, operates between checkpoints and then returns
-- .............................................................................
local function run_simulation(pinit, configure_mara, runargs)

   local function InitSimulation(pinit, configure_mara)
      configure_mara()
      init_prim(pinit)
      boundary.ApplyBoundaries()

      local Status = { }

      Status.CurrentTime = 0.0
      Status.Iteration   = 0
      Status.Checkpoint  = 0
      Status.Timestep    = 0.0

      local datadir = string.format("data/%s", runargs.id)
      os.execute(string.format("mkdir -p %s", datadir))
      os.execute(host.Filesystem(datadir))

      if not runargs.quiet then print_mara() end
      return Status
   end

   local function HandleErrors(Status, attempt)
      return 0
   end

   local Status = InitSimulation(pinit, configure_mara)

   local t0 = Status.CurrentTime
   local attempt = 0

   if runargs.dim == 2 and not runargs.noplot then
      visual.open_window()
   end

   while Status.CurrentTime - t0 < runargs.tmax do

      if runargs.dim == 2 and not runargs.noplot then
         local prim = get_prim()
         if runargs.dim == 2 then
            local draw_array = prim.rho
            visual.draw_texture(draw_array)
         else
            local draw_array = prim.rho[':,:,'..runargs.N/2]
            visual.draw_texture(draw_array)
         end
      end

      if HandleErrors(Status, attempt) ~= 0 then
         return 1
      end
      attempt = attempt + 1

      local dt = Status.Timestep
      local kzps, errors = advance(dt)

      if errors == 0 then
         driving.Advance(dt)
         driving.Resample()

         if not runargs.quiet then
            print(string.format("%05d(%d): t=%5.4f dt=%5.4e %3.2fkz/s %3.2fus/(z*Nq)",
                                Status.Iteration, attempt-1, Status.CurrentTime, dt,
                                kzps, 1e6/8/(1e3*kzps)))
            io.flush()
         end

         attempt = 0
         Status.CurrentTime = Status.CurrentTime + Status.Timestep
         Status.Timestep = get_timestep(runargs.CFL)
         Status.Iteration = Status.Iteration + 1
      end
   end
   return Status
end


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

-- -----------------------------------------------------------------------------
-- http://lua-users.org/wiki/SplitJoin
-- -----------------------------------------------------------------------------
function string_split(self, sSeparator, nMax, bRegexp)
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


local function parse_args(runargs)
   for k,v in pairs(cmdline.opts) do
      if type(runargs[k]) == 'number' then
         runargs[k] = tonumber(v)
      else
         runargs[k] = v
      end
   end
   print("command line arguments:")
   for k,v in pairs(runargs) do
      print(k,v)
   end
end

return { deepcopy=deepcopy,
         plot=plot,
         parse_args=parse_args,
         run_simulation=run_simulation,
	 string_split=string_split }
