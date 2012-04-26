
local function do_units()
   local LIGHT_SPEED = 2.99792458000e+10 -- cm/s

   local Density = 1e13         -- gm/cm^3
   local V       = LIGHT_SPEED  -- cm/s
   local Length  = 1e2          -- cm
   local Mass    = Density * Length^3.0
   local Time    = Length / V
   set_units(Length, Mass, Time)
   units.Print()
end

local function LoadMicroPh(fname)
   print("Loading tabulated EOS from " .. fname)
   print("[note] be sure you have already set the physics units before calling")

   h5_open_file(fname, "r")
   local eos_terms = { }

   for _,v in pairs({"pressure", "internal_energy",
		     "sound_speed", "density", "temperature"}) do
      eos_terms[v] = h5_read_array(v)
   end
   h5_close_file()

   local D = eos_terms["density"    ][':,0'] * units.GramsPerCubicCentimeter()
   local T = eos_terms["temperature"]['0,:']

   local p = eos_terms["pressure"] * units.MeVPerCubicFemtometer()
   local u = eos_terms["internal_energy"] * units.MeVPerCubicFemtometer()
   local c = eos_terms["sound_speed"]

   set_eos("tabulated", {D=D, T=T, p=p, u=u, c=c})
end

local _tabeos_ = { LoadMicroPh=LoadMicroPh,
		   MakeNeutronStarUnits=do_units }
return _tabeos_
