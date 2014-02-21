
if arg[2] == '--explain' then
   print "3d driven turbulence research problem"
   return
end

local MPI     = require 'MPI'
local cow     = require 'cow'
local hdf5    = require 'lua-hdf5.LuaHDF5'
local Mara    = require 'Mara'
local unigrid = require 'new-unigrid'
local array   = require 'array'
local json    = require 'json'
local util    = require 'util'

local RunArgs = {
   N       = 16,
   id      = "test",
   cool    = "E4",     -- for (internal energy)^4 or T4, for temperature^4, or none
   coolP   = 0.1,      -- cooling parameter, reference temp(T4) or e(E4)
   zeta    = 1.0,      -- driving vorticity (0=none -> 1=full)
   F0      = 0.01,     -- driving field power coefficient
   B0      = 1e10,     -- B-field, in Gauss
   D0      = 1.00,     -- density, in code units
   P0      = 0.05,     -- pressure, in code units
   CFL     = 0.24,
   tmax    = 72.0,
   cpi     = 1.0,      -- checkpoint interval, in code time units
   restart = "none",
   rmdivB  = false,    -- run a div-clean on the input B-field
   eosfile = "none",   -- i.e. nseos.h5
   fluid   = "rmhd",   -- euler, srhd, rmhd
   gamma   = 4/3,      -- adiabatic gamma (1.001 for isothermal)
   pspec   = 0,        -- pspec > 0 take power spectra every that many iterations
   pdfs    = 0,        -- pdfs > 0 take PDF's ' '
   problem = "drvtrb", -- or KH
   drive   = true,     -- set to false to disable driving,
   scheme  = 'weno',   -- weno or hllc,
}

local PowerSpectrumFile = ""
local DistributionsFile = ""



-- *****************************************************************************
-- Process command line options
-- .............................................................................
local function process_cmdline()
   for k,v in pairs(arg) do
      local eqpos = string.find(v, '=')
      if eqpos then
         local key = string.sub(v, 0, eqpos-1)
         local val = string.sub(v, eqpos+1)

         if type(RunArgs[key]) == 'boolean' then
            RunArgs[key] = (val == "1" or val == "true")
         elseif type(RunArgs[key]) == 'number' then
            RunArgs[key] = tonumber(val)
         elseif type(RunArgs[key]) == 'string' then
            RunArgs[key] = val
         else
            print("[Mara] warning: ignoring unrecognized option "..key)
         end
      end
   end
end




-- *****************************************************************************
-- Take one-point pdf's
-- .............................................................................
local function Distributions(primitive, gname)
   local fname = DistributionsFile

   if primitive:domain():get_rank() == 0 then
      if not util.file_exists(fname) then
         local h5f = hdf5.File(fname, 'w')
         h5f:close()
      end
   end

   local P = primitive._dfield
   cow.srhdpack.onepointpdfs(P, 'proper-rho', fname, gname, 1e-6, 1e6)
   cow.srhdpack.onepointpdfs(P, 'gamma-rho' , fname, gname, 1e-6, 1e6)
   cow.srhdpack.onepointpdfs(P, 'gamma-beta', fname, gname, 1e-6, 1e6)
end




-- *****************************************************************************
-- Project out divergence of B
-- .............................................................................
local function RemoveDivergenceB(primitive)
   local Bgrid = unigrid.UnigridDataField(primitive:domain(), {'Bx','By','Bz'})
   local P = primitive:array()
   local B = Bgrid:array()
   B[{nil,nil,nil,{0,1}}] = P[{nil,nil,nil,{5,6}}]
   B[{nil,nil,nil,{1,2}}] = P[{nil,nil,nil,{6,7}}]
   B[{nil,nil,nil,{2,3}}] = P[{nil,nil,nil,{7,8}}]
   Bgrid:project_out_divergence()
   P[{nil,nil,nil,{5,6}}] = B[{nil,nil,nil,{0,1}}]
   P[{nil,nil,nil,{6,7}}] = B[{nil,nil,nil,{1,2}}]
   P[{nil,nil,nil,{7,8}}] = B[{nil,nil,nil,{2,3}}]
end




-- *****************************************************************************
-- Take power spectrum
-- .............................................................................
local function PowerSpectrum(primitive, which, gname, helmholtz)
   local start = os.clock()
   local binloc, binval
   local fname = PowerSpectrumFile
   if which == 'velocity' then
      local velocity = unigrid.UnigridDataField(primitive:domain(), {'vx','vy','vz'})
      local P = primitive:array()
      local V = velocity:array()
      V[{nil,nil,nil,{0,1}}] = P[{nil,nil,nil,{2,3}}]
      V[{nil,nil,nil,{1,2}}] = P[{nil,nil,nil,{3,4}}]
      V[{nil,nil,nil,{2,3}}] = P[{nil,nil,nil,{4,5}}]
      if helmholtz == 'solenoidal' then
	 velocity:project_out_divergence()
	 which = which .. "-solenoidal"
      elseif helmholtz == 'dilatational' then
	 velocity:project_out_curl()
	 which = which .. "-dilatational"
      end
      binloc, binval = velocity:power_spectrum(128)
   elseif which == 'magnetic' then
      local magnetic = unigrid.UnigridDataField(primitive:domain(), {'Bx','By','Bz'})
      local P = primitive:array()
      local B = magnetic:array()
      B[{nil,nil,nil,{0,1}}] = P[{nil,nil,nil,{5,6}}]
      B[{nil,nil,nil,{1,2}}] = P[{nil,nil,nil,{6,7}}]
      B[{nil,nil,nil,{2,3}}] = P[{nil,nil,nil,{7,8}}]
      if helmholtz == 'solenoidal' then
	 magnetic:project_out_divergence()
	 which = which .. "-solenoidal"
      elseif helmholtz == 'dilatational' then
	 magnetic:project_out_curl()
	 which = which .. "-dilatational"
      end
      binloc, binval = magnetic:power_spectrum(128)
   elseif which == 'kinetic' then
      local kinetic = unigrid.UnigridDataField(primitive:domain(), {'Kx','Ky','Kz'})
      local P = primitive:array()
      local K = kinetic:array()
      local S = K:shape()
      local vx = array.array{S[1], S[2], S[3]}
      local vy = array.array{S[1], S[2], S[3]}
      local vz = array.array{S[1], S[2], S[3]}
      local D0 = array.array{S[1], S[2], S[3]}
      D0[nil] = P[{nil,nil,nil,{0,1}}]
      vx[nil] = P[{nil,nil,nil,{2,3}}]
      vy[nil] = P[{nil,nil,nil,{3,4}}]
      vz[nil] = P[{nil,nil,nil,{4,5}}]
      local D0_vec = D0:vector()
      local vx_vec = vx:vector()
      local vy_vec = vy:vector()
      local vz_vec = vz:vector()
      for i=0,#D0_vec-1 do
         local a = (0.5 * D0_vec[i])^0.5
         vx_vec[i] = vx_vec[i] * a
         vy_vec[i] = vy_vec[i] * a
         vz_vec[i] = vz_vec[i] * a
      end
      K[{nil,nil,nil,{0,1}}] = vx
      K[{nil,nil,nil,{1,2}}] = vy
      K[{nil,nil,nil,{2,3}}] = vz
      binloc, binval = kinetic:power_spectrum(128)
   end
   if primitive:domain():get_rank() == 0 then
      local h5f = hdf5.File(fname, 'r+')
      local grp = hdf5.Group(h5f, gname)
      local sgr = hdf5.Group(grp, which)
      sgr['binloc'] = binloc
      sgr['binval'] = binval
      h5f:close()
   end
   print(string.format("[Mara] PowerSpectrum took %f seconds",
                       os.clock() - start))
end


-- *****************************************************************************
-- Function to collect all available measurements from Mara
-- .............................................................................
local function Measurements(Status)
   local meas = { }
   local dtmeas = Status.CurrentTime - Status.LastMeasurementTime

   meas.U                    = Mara.measure.mean_cons()
   meas.P                    = Mara.measure.mean_prim()

   meas.energies             = Mara.measure.mean_energies()
   meas.mean_velocity        = Mara.measure.mean_velocity()

   meas.mean_T, meas.max_T   = Mara.measure.mean_max_temperature()
   meas.mean_B, meas.max_B   = Mara.measure.mean_max_magnetic_field()
   meas.mean_Ms, meas.max_Ms = Mara.measure.mean_max_sonic_mach()
   meas.mean_Ma, meas.min_Ma = Mara.measure.mean_min_alfvenic_mach()
   meas.max_lorentz_factor   = Mara.measure.max_lorentz_factor()
   meas.cooling_rate         = Mara.cooling_rate(dtmeas)

   meas.Status = util.deepcopy(Status)
   return meas
end


-- *****************************************************************************
-- Error handlers for different fluids
-- .............................................................................
local function HandleErrorsEuler(P, Status, attempt)
   Mara.set_advance("rk4")
   if attempt == 0 then -- healthy time-step
      set_riemann("hllc")
      set_godunov("plm-split")
      Status.Timestep = 1.0 * Status.Timestep
      return 0
   else
      return 1
   end
end

local function HandleErrorsSrhd(P, Status, attempt)
   if attempt == 0 then -- healthy time-step
      Mara.set_advance("rk3")
      if RunArgs.scheme == 'weno' then
         Mara.set_godunov("weno-split")
      elseif RunArgs.scheme == 'hllc' then
         Mara.set_godunov("plm-split")
         Mara.set_riemann("hllc")
      else
         error('no such scheme: '..RunArgs.scheme)
      end
      Status.Timestep = 1.0 * Status.Timestep
      return 0
   elseif attempt == 1 then
      Status.Timestep = 0.5 * Status.Timestep
      Mara.diffuse(P, 0.2)
      return 0
   elseif attempt == 2 then
      Status.Timestep = 0.5 * Status.Timestep
      Mara.diffuse(P, 0.2)
      return 0
   elseif attempt == 3 then
      Mara.set_godunov("plm-split")
      Mara.set_riemann("hll")
      Status.Timestep = 0.5 * Status.Timestep
      Mara.diffuse(P, 0.2)
      return 0
   elseif attempt == 4 then
      Status.Timestep = 0.5 * Status.Timestep
      Mara.diffuse(P, 0.2)
      return 0
   else
      return 1
   end
end

local function HandleErrorsRmhd(P, Status, attempt)
   Mara.set_advance("single")
   if attempt == 0 then -- healthy time-step
      Mara.set_godunov("plm-muscl")
      Mara.set_riemann("hlld")
      Mara.config_solver({theta=2.0}, true)
      Status.Timestep = 1.0 * Status.Timestep
      return 0
   elseif attempt == 1 then
      Mara.set_godunov("plm-muscl")
      Mara.config_solver({theta=1.25}, true)
      Mara.diffuse(P, 0.2)
      Status.Timestep = 0.5 * Status.Timestep
      return 0
   elseif attempt == 2 then
      Mara.set_godunov("plm-muscl")
      Mara.config_solver({theta=1.0}, true)
      Mara.diffuse(P, 0.2)
      Status.Timestep = 0.5 * Status.Timestep
      return 0
   elseif attempt == 3 then
      Mara.set_godunov("plm-muscl")
      Mara.config_solver({theta=0.0}, true)
      Mara.set_riemann("hll")
      Mara.diffuse(P, 0.2)
      Status.Timestep = 0.5 * Status.Timestep
      return 0
   elseif attempt == 4 then
      Mara.diffuse(P, 0.2)
      Status.Timestep = 0.5 * Status.Timestep
      return 0
   elseif attempt == 5 then
      Mara.diffuse(P, 0.2)
      Status.Timestep = 0.5 * Status.Timestep
      return 0
   else
      return 1
   end
end

local function HandleErrorsCascade(P, Status, attempt)
   Mara.set_advance("rk3")
   Mara.set_riemann("hll")
   Mara.set_godunov("plm-split")
   Mara.config_solver({theta=1.5}, true)
   if attempt == 0 then -- healthy time-step
      return 0
   elseif attempt == 1 then
      Status.Timestep = 0.5 * Status.Timestep
      Mara.diffuse(P, 0.5)
      return 0
   elseif attempt == 2 then
      Status.Timestep = 0.5 * Status.Timestep
      Mara.diffuse(P, 0.5)
      return 0
   elseif attempt == 3 then
      Status.Timestep = 0.5 * Status.Timestep
      Mara.diffuse(P, 0.5)
      return 0
   else
      return 1
   end
end


-- *****************************************************************************
-- Units for neutron star problem
-- .............................................................................
local function MakeNeutronStarUnits()
   local LIGHT_SPEED = 2.99792458000e+10 -- cm/s
   local Density = 1e13         -- gm/cm^3
   local V       = LIGHT_SPEED  -- cm/s
   local Length  = 1e2          -- cm
   local Mass    = Density * Length^3.0
   local Time    = Length / V
   Mara.set_units(Length, Mass, Time)
end


-- *****************************************************************************
-- Load tabulated equation of state
-- .............................................................................
local function LoadMicroPh(fname)
   print("[Mara] loading tabulated EOS from "..fname)
   print("[Mara] note: be sure you have already set the physics units")

   local eos_file = hdf5.File(fname, 'r')

   -- Load EOS data from the file
   local D = eos_file["density"    ][{{nil},{0,1}}]:table()
   local T = eos_file["temperature"][{{0,1},{nil}}]:table()
   local p = eos_file["pressure"       ][nil]:table()
   local u = eos_file["internal_energy"][nil]:table()
   local c = eos_file["sound_speed"    ][nil]:table()

   -- Convert EOS from physics to code units (EOS is stored in NS units)
   for i=1,#D do D[i] = D[i] * Mara.units.GramsPerCubicCentimeter() end
   for i=1,#p do p[i] = p[i] * Mara.units.MeVPerCubicFemtometer() end
   for i=1,#u do u[i] = u[i] * Mara.units.MeVPerCubicFemtometer() end

   Mara.set_eos("tabulated", {D=D, T=T, p=p, u=u, c=c})
   eos_file:close()
end


-- *****************************************************************************
-- Main driver, operates between checkpoints and then returns
-- .............................................................................
local function RunSimulation(Primitive, Status, MeasureLog, Howlong,
			     HandleErrors)

   local t0 = Status.CurrentTime
   local attempt = 0
   local NumberOfConserved = { euler=5, srhd=5, rmhd=8 }

   while Status.CurrentTime - t0 < Howlong do
      collectgarbage()

      local P = Primitive:array():buffer()
      local stopfname = string.format("data/%s/MARA_STOP", RunArgs.id)
      local stopfile = io.open(stopfname, "r")
      if stopfile then
         print("Mara exiting, detected", stopfname)
         stopfile.close()
         mpi_barrier()
         if mpi_get_rank() == 0 then
            os.remove(stopfname)
         end
         return "STOP"
      end

      -- Measurements are made at the beginning of the timestep
      -- .......................................................................
      if attempt == 0 then
         if Status.Iteration % 10 == 0 then
            Mara.set_primitive(P)
            MeasureLog[Status.Iteration] = Measurements(Status)
            Status.LastMeasurementTime = Status.CurrentTime
         end

         if Status.Iteration % RunArgs.pspec == 0 and RunArgs.pspec ~= 0 then
            local gname = string.format("pspec-%05d", Status.Iteration)
            PowerSpectrum(Primitive, 'magnetic', gname, 'solenoidal')
            PowerSpectrum(Primitive, 'magnetic', gname, 'dilatational')
            PowerSpectrum(Primitive, 'velocity', gname, 'solenoidal')
            PowerSpectrum(Primitive, 'velocity', gname, 'dilatational')
         end

         if Status.Iteration % RunArgs.pdfs == 0 and RunArgs.pdfs ~= 0 then
            local gname = string.format("pdf-%05d", Status.Iteration)
            Distributions(Primitive, gname)
         end
      end

      -- 'attempt' == 0 when the previous iteration completed without errors
      -- .......................................................................
      if HandleErrors(P, Status, attempt) ~= 0 then
         return 1
      end
      attempt = attempt + 1

      local dt = Status.Timestep
      local kzps, errors = Mara.advance(P, dt)

      if errors == 0 then
         if RunArgs.problem == "drvtrb" and RunArgs.drive then
            Mara.driving.Advance(dt)
            if Status.Iteration % 10 == 0 then
               Mara.driving.Resample()
            end
         end

         local Nq = NumberOfConserved[RunArgs.fluid]
         print(string.format("%05d(%d): t=%5.4f dt=%5.4e %3.1fkz/s %3.2fus/(z*Nq)",
                             Status.Iteration, attempt-1, Status.CurrentTime, dt,
                             kzps, 1e6/Nq/(1e3*kzps)))
         io.flush()

         attempt = 0
         Status.Timestep = Mara.get_timestep(RunArgs.CFL)
         Status.CurrentTime = Status.CurrentTime + Status.Timestep
         Status.Iteration = Status.Iteration + 1
      end
   end
   return 0
end


-- *****************************************************************************
-- Write restart files (checkpoints)
-- .............................................................................
local function CheckpointWrite(Primitive, Status, MeasureLog, OptionalName)
   local datadir = string.format("data/%s", RunArgs.id)
   local chkpt
   if OptionalName then
      chkpt = string.format("%s/chkpt.%s.h5", datadir, OptionalName)
   else
      Status.Checkpoint = Status.Checkpoint + 1
      chkpt = string.format("%s/chkpt.%04d.h5", datadir, Status.Checkpoint)
   end

   if Primitive:domain():get_rank() == 0 then
      if util.file_exists(chkpt) then
         os.execute(string.format('rm %s', chkpt))
      end
   end

   Primitive:write(chkpt, 'prim')

   if Primitive:domain():get_rank() == 0 then
      local version = Mara.git_sha()
      local program = " "
      local f = io.open(string.sub(debug.getinfo(1).source, 2), 'r')
      if f then
         program = f:read('*all')
         f:close()
      end
      local chkpt_h5 = hdf5.File(chkpt, 'r+')
      chkpt_h5["runargs"] = json.encode(RunArgs)
      chkpt_h5["measure"] = json.encode(MeasureLog)
      chkpt_h5["program"] = program
      chkpt_h5["version"] = version
      if RunArgs.problem == "drvtrb" and RunArgs.drive then
         chkpt_h5["driving"] = json.encode(Mara.driving.Serialize())
      end
      local status = hdf5.Group(chkpt_h5, 'status')
      for k,v in pairs(Status) do
         status[k] = v
      end
      status:close()
      chkpt_h5:close()
   end
end


-- *****************************************************************************
-- Read restart files (checkpoints)
-- .............................................................................
local function CheckpointRead(chkpt, Primitive)

   Primitive:read(chkpt, 'prim')

   local chkpt_h5 = hdf5.File(chkpt, 'r')
   local status_h5 = hdf5.Group(chkpt_h5, "status")
   local MeasureLog = json.decode(chkpt_h5["measure"]:value())

   -- workaround! json-ing a table which was previously empty doesn't work?
   if chkpt_h5["measure"]:value() == '[]' then
      print('resetting measure to { }')
      MeasureLog = { }
   end

   local Status = { }
   for key,val in pairs(status_h5) do
      Status[key] = val:value()
   end
   status_h5:close()
   chkpt_h5:close()
   return Status, MeasureLog
end


-- *****************************************************************************
-- Initialize and run simulation
-- .............................................................................
local function main()
   process_cmdline()

   MPI.Init()
   cow.init(0, nil, 0) -- to reopen stdout to dev/null
   Mara.start()
   Mara.set_fluid(RunArgs.fluid)
   Mara.set_advance('single')
   Mara.set_godunov('plm-muscl')
   Mara.set_boundary('periodic')
   Mara.set_riemann('hlld')
   Mara.config_solver{IS="sz10", sz10A=100.0}

   MakeNeutronStarUnits() -- must come before loading the EOS
   if RunArgs.eosfile == "none" then
      Mara.set_eos('gamma-law', RunArgs.gamma)
   else
      LoadMicroPh('nseos.h5')
   end

   local prim_names = Mara.fluid.GetPrimNames()
   local Nq = #prim_names
   local Ng = 3
   local Nx, Ny, Nz, L0, L1, pinit, HandleErrors

   if RunArgs.problem == "drvtrb" then
      Nx = RunArgs.N
      Ny = RunArgs.N
      Nz = RunArgs.N
      L0 = { -0.5, -0.5, -0.5 }
      L1 = {  0.5,  0.5,  0.5 }

      if RunArgs.cool == "none" then
         print("[drvtrb] disabling cooling")
      elseif RunArgs.cool == "T4" then
         Mara.set_cooling("T4", 25.0, 100.0)
      elseif RunArgs.cool == "E4" then
         Mara.set_cooling("E4", RunArgs.coolP, 100.0)
      end

      pinit = function(x,y,z)
         local B0 = RunArgs.B0 * Mara.units.Gauss()
         local D0 = RunArgs.D0
         local P0 = RunArgs.P0
         return { D0, P0, 0, 0, 0, B0, 0.0, 0.0 }
      end
      HandleErrors = ({ euler=HandleErrorsEuler,
			srhd=HandleErrorsSrhd,
			rmhd=HandleErrorsRmhd })[RunArgs.fluid]

   elseif RunArgs.problem == "KH" then
      Nx = RunArgs.N
      Ny = RunArgs.N * 2
      Nz = RunArgs.N
      L0 = { -0.5, -1.0, -0.5 }
      L1 = {  0.5,  1.0,  0.5 }

      pinit = function(x,y,z)
         local Vs = 0.5 -- Shearing velocity
         local a = 0.01 -- Width of shearing layer
         local A0 = 0.1 -- Sinusoidal perturbation amplitude
         local sig = 0.1 -- Length scale of its decay
         local sin, tanh, exp, pi = math.sin, math.tanh, math.exp, math.pi
         local B0 = 1e-3
         local D0 = math.abs(y) < 0.5 and 0.01 or 1.0
         local P0 = 1.0
         local vx = y > 0.0 and Vs*tanh((y-0.5)/a) or -Vs*tanh((y+0.5)/a)
         local vy = (y > 0.0 and
                     A0*Vs*sin(2*pi*x)*exp(-((y-0.5)/sig)^2) or
                     A0*Vs*sin(2*pi*x)*exp(-((y+0.5)/sig)^2)*-1)
         local vz = (math.random() - 0.5) * 0.01
         return { D0, P0, vx, vy, vz, B0, 0.0, 0.0 }
      end

      HandleErrors = ({ euler=HandleErrorsEuler,
			srhd=HandleErrorsSrhd,
			rmhd=HandleErrorsRmhd })[RunArgs.fluid]

   elseif RunArgs.problem == "cascade" then
      Nx = RunArgs.N
      Ny = RunArgs.N
      Nz = RunArgs.N
      L0 = { -0.5, -0.5, -0.5 }
      L1 = {  0.5,  0.5,  0.5 }
      pinit = function(x,y,z)
         local D0 = RunArgs.D0
         local P0 = RunArgs.P0
         return { D0, P0, 0, 0, 0, 0.0, 0.0, 0.0 }
      end
      HandleErrors = HandleErrorsCascade
   else
      error("[drvtrb] problem must be either drvtrb, KH, or cascade")
   end

   local domain = unigrid.UnigridDomain({Nx, Ny, Nz}, Ng)
   local primitive = unigrid.UnigridDataField(domain, prim_names)
   local domain_comm = domain:get_comm()
   local P = primitive:array()

   if hdf5.have_mpio() then
      domain:set_collective(1)
   end

   math.randomseed(domain:get_rank())
   Mara.set_domain(L0, L1, {Nx, Ny, Nz}, Nq, Ng, domain_comm)

   local MeasureLog = { }
   local Status = { }

   if RunArgs.restart == "none" then
      Mara.init_prim(P:buffer(), pinit) -- start a new model from scratch
      Status.CurrentTime = 0.0
      Status.Iteration   = 0
      Status.Checkpoint  = 0
      Status.Timestep    = 0.0
      Status.LastMeasurementTime = 0.0

      if RunArgs.problem == "drvtrb" and RunArgs.drive then
         print("[drvtrb] enabling driving field")
         local field = Mara.new_ou_field(3, RunArgs.F0, RunArgs.zeta, 3, 12345)
         Mara.set_driving(field)
      else
         print("[drvtrb] disabling driving field")
      end
   else
      if RunArgs.problem == "drvtrb" and RunArgs.drive then
         print("[drvtrb] enabling serialized driving field")
         local restart_file = hdf5.File(RunArgs.restart, 'r')
         local field = json.decode(restart_file["driving"]:value())
         Mara.set_driving(field)
         restart_file:close()
      end
      Status, MeasureLog = CheckpointRead(RunArgs.restart, primitive)
      if RunArgs['rmdivB'] then
	 RemoveDivergenceB(primitive)
      end
   end

   PowerSpectrumFile = string.format("data/%s/%s.spec.h5", RunArgs.id, RunArgs.id)
   DistributionsFile = string.format("data/%s/%s.pdfs.h5", RunArgs.id, RunArgs.id)
   local datadir = string.format("data/%s", RunArgs.id)

   if domain:get_rank() == 0 then
      os.execute(string.format("mkdir -p %s", datadir))
      if RunArgs.pspec ~= 0 then
         local testf = io.open(PowerSpectrumFile, "r")
         if not testf then
            local h5f = hdf5.File(PowerSpectrumFile, "w")
            h5f:close()
         else
            testf:close()
         end
      end
   end

   Mara.show()
   Mara.units.Print()
   print('\n\t***************************')
   print('\tRuntime arguments:')
   print('\t***************************')
   util.pretty_print(RunArgs, '\t')
   print('\t***************************\n')


   while Status.CurrentTime < RunArgs.tmax do
      local error = RunSimulation(primitive, Status, MeasureLog, RunArgs.cpi,
				  HandleErrors)
      if error == "STOP" then
         print("exiting upon request\n")
         CheckpointWrite(primitive, Status, MeasureLog, "stop")
         break
      elseif error ~= 0 then
         print("exiting due to failures\n")
         CheckpointWrite(primitive, Status, MeasureLog, "fail")
         break
      end
      CheckpointWrite(primitive, Status, MeasureLog)
   end

   primitive = nil
   domain = nil
   collectgarbage()

   Mara.close()
   MPI.Finalize()
end

main()
