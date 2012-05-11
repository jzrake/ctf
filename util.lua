

local host = require 'host'
local json = require 'json'
local util = { }

-- *****************************************************************************
-- Main driver, operates between checkpoints and then returns
-- .............................................................................
function util.run_simulation(pinit, configure_mara, runargs)

   local function InitSimulation(pinit, configure_mara)
      configure_mara()
      init_prim(pinit)
      boundary.ApplyBoundaries()

      local Status = { }

      Status.CurrentTime = 0.0
      Status.Iteration   = 0
      Status.Checkpoint  = 0
      Status.LastCheckpointTime = 0.0
      Status.Timestep    = 0.0

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

      if (runargs.cpi > 0.0 and
	  Status.CurrentTime - Status.LastCheckpointTime > runargs.cpi) then
	 util.write_checkpoint(Status, runargs)
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
         Status.Iteration = Status.Iteration + 1
	 if runargs.fixdt < 0.0 then
	    Status.Timestep = get_timestep(runargs.CFL)
	 else
	    Status.Timestep = runargs.fixdt
	 end
      end
   end
   return Status
end


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
-- Read and write checkpoints with HDF5
-- .............................................................................
function util.write_checkpoint(Status, RunArgs)
   Status.Checkpoint = Status.Checkpoint + 1
   Status.LastCheckpointTime = Status.CurrentTime

   local datadir = string.format("data/%s", RunArgs.id)
   local version = mara_version()
   local chkpt = string.format("%s/chkpt.%04d.h5", datadir, Status.Checkpoint)

   if mpi_get_rank() == 0 then
      os.execute(string.format("mkdir -p %s", datadir))
      os.execute(host.Filesystem(datadir))
      h5_open_file(chkpt, "w")
      h5_write_numeric_table("status", Status)
      h5_write_string("runargs", json.encode(RunArgs))
      h5_write_string("version", version)
      h5_close_file()
   end
   write_prim(chkpt, host.CheckpointOptions)
end

function util.read_checkpoint(chkpt)
   h5_open_file(chkpt, "r")
   local status = h5_read_numeric_table("status")
   h5_close_file()
   read_prim(chkpt, host.CheckpointOptions)
   return status
end

-- *****************************************************************************
-- Function to call Gnuplot from Lua using popen
-- .............................................................................
function util.plot(series, tpause)
   local gp = io.popen("gnuplot", 'w')

   if util.RunArgs.pdf then
      gp:write("set terminal postscript enhanced color\n")
      gp:write(string.format("set output '| ps2pdf - %s.pdf'\n", util.RunArgs.id))
      gp:write(string.format("set title '%s'\n", util.RunArgs.id))
   end

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
   if not util.RunArgs.pdf then
      gp:write(string.format("pause %f\n", tpause or 100.0))
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

-- -----------------------------------------------------------------------------
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
