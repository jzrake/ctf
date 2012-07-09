

local tests = { }
tests.RunArgs = { }

tests.IsentropicPulse = {
   get_pinit =
      function(self)
      local g = function(x,y,z)
         local L = 1.0
         local n = self.mode
         local K = self.entropy_ref
         local Gamma = self.Gamma
         local rho_ref = self.rho_ref
         local pre_ref = K * rho_ref ^ Gamma
         local cs_ref = (Gamma * pre_ref / rho_ref)^0.5

         local function f(x)
            return math.sin(n*math.pi*x/L)^2
         end

         local rho = rho_ref * (1.0 + f(x))
         local pre = K * rho ^ Gamma
         local cs = (Gamma * pre/rho)^0.5
         local vx = 2 / (Gamma - 1) * (cs - cs_ref)
         return { rho, pre, vx, 0, 0 }
      end
      return g
      end,
      entropy =
         function(self, rho, pre)
         return pre / rho^self.Gamma
         end,
         entropy_ref = 0.1, -- set equal to K above
         mode = 2,
         Gamma = 1.4,
         rho_ref = 1.0
}

tests.Explosion = {
   get_pinit =
      function(self, time)
      local g =
         function(x,y,z)
         local dim = tests.RunArgs.dim
         if dim >= 1 then x = x - 0.5 end
         if dim >= 2 then y = y - 0.5 end
         if dim >= 3 then z = z - 0.5 end
         local r2 = x*x + y*y + z*z
         if r2 < 0.01 then
            return { 1.000, 1.0, 0, 0, 0, self.Bx, 0, 0 }
         else
            return { 0.125, 0.1, 0, 0, 0, self.Bx, 0, 0 }
         end
         end
         return g
      end,
      Bx = 4.0
}

tests.KelvinHelmholtz = {
   get_pinit = function(self)
      local g =
         function(x,y,z)
         local x = x - 0.5
         local y = y - 0.5
         local z = z - 0.5
	 local vx = 0.02*(math.random() - 0.5) + self:vx_profile(x,y,z)
         local vy = 0.02*(math.random() - 0.5)
	 local rho = self:rho_profile(x,y,z)
         return { rho, 2.5, vx, vy, 0.0 }
         end
         return g
   end,
   rho_profile = function(self,x,y,z)
      return math.abs(y) > 0.25 and 1.0 or 2.0
   end,
   vx_profile = function(self,x,y,z)
      return math.abs(y) > 0.25 and -.5 or 0.5
   end
}


tests.SmoothKelvinHelmholtz = {
   get_pinit =
      function(self)
      local g =
         function(x,y,z)
         local P0 = 2.5
         local rho1 = 1.0
         local rho2 = 2.0
         local L = .025
         local U1 = 0.5
         local U2 = -0.5
         local w0 = 0.01
         local sin, exp = math.sin, math.exp

         local vy = w0*sin(4*math.pi*x)
         local rho,vx
         if y < 0.25 then
            rho = rho1 - 0.5*(rho1-rho2)*exp( (y-0.25)/L)
            vx  = U1   - 0.5*( U1 - U2 )*exp( (y-0.25)/L)
         elseif y < 0.5 then
            rho = rho2 + 0.5*(rho1-rho2)*exp(-(y-0.25)/L)
            vx  = U2   + 0.5*( U1 - U2 )*exp(-(y-0.25)/L)
         elseif y < 0.75 then
            rho = rho2 + 0.5*(rho1-rho2)*exp( (y-0.75)/L)
            vx  = U2   + 0.5*( U1 - U2 )*exp( (y-0.75)/L)
         else
            rho = rho1 - 0.5*(rho1-rho2)*exp(-(y-0.75)/L)
            vx  = U1   - 0.5*( U1 - U2 )*exp(-(y-0.75)/L)
         end
         return { rho, P0, vx, vy, 0.0 }
         end
         return g
      end
}

tests.DensityWave = {
   get_pinit =
      function(self, time)
      local g =
         function(x,y,z)
         return {self:true_solution(time or 0.0, x, y, z),
                 1.0,
                 self.velocity[1],
                 self.velocity[2],
                 self.velocity[3]}
         end
         return g
      end,
      true_solution =
         function(self, t, x, y, z)
         local k = self.mode
         local v = self.velocity
         local kdotx = (x - v[1]*t)*k[1] + (y - v[2]*t)*k[2] + (z - v[3]*t)*k[3]
         return self.rho_ref*(1.0 + self.eps*math.cos(2*math.pi*(kdotx)))
         end,
         velocity = { 1.0, 0.0, 0.0 },
         mode = { 1, 0, 0 },
         eps = 3.2e-1,
         rho_ref = 1.0
}

tests.CollidingShocks = {
   get_pinit =
      function(self, time)
      local g =
         function(x,y,z)
         if x < 0.1 then
            return { 1.0, 1e+3, 0, 0, 0 }
         elseif 0.1 < x and x < 0.9 then
            return { 1.0, 1e-2, 0, 0, 0 }
         else
            return { 1.0, 1e+2, 0, 0, 0 }
         end
         end
         return g
      end
}

tests.TangentialVelocity = {
   get_pinit =
      function(self, time)
      local g =
         function(x,y,z)
         local x = x - 0.5
         local vy = math.exp(-x*x / 0.01)
         return { 1.0, 1.0, 0.0, vy, 0.0 }
         end
         return g
      end
}


local make_2state_setup =
   function(states, opts)
   return {
      get_pinit =
         function(self)
         local g =
            function(x,y,z)

            local n = load("return"..tests.RunArgs.angle)()
            local nmag = (n[1]*n[1] + n[2]*n[2] + n[3]*n[3])^0.5
            local xn = (n[1]*x + n[2]*y + n[3]*z) / nmag

            if not (opts and opts.reverse) then
               if xn < 0.5 then
                  return states.Pl
               else
                  return states.Pr
               end
            else
               if xn < 0.5 then
                  return states.Pr
               else
                  return states.Pl
               end
            end
            end
            return g
         end
          }
   end

   tests.SrhdCase1_DFIM98 = make_2state_setup
   ({ Pl = { 10.0, 13.30, 0.0, 0.0, 0.0 },
      Pr = {  1.0,  1e-6, 0.0, 0.0, 0.0 } })

   tests.SrhdCase2_DFIM98 = make_2state_setup
   ({ Pl = { 1, 1e+3, 0.0, 0.0, 0.0 },
      Pr = { 1, 1e-2, 0.0, 0.0, 0.0 } })

   tests.SrhdHardTransverse_RAM = make_2state_setup
   ({ Pl = { 1, 1e+3, 0.0, 0.9, 0.0 },
      Pr = { 1, 1e-2, 0.0, 0.9, 0.0 } })

   tests.Shocktube1 = make_2state_setup
   ({ Pl = { 1.000, 1.000, 0.000, 0.0, 0.0 },
      Pr = { 0.125, 0.100, 0.000, 0.0, 0.0 } })

   tests.Shocktube2 = make_2state_setup
   ({ Pl = { 1.000, 0.400,-2.000, 0.0, 0.0 },
      Pr = { 1.000, 0.400, 2.000, 0.0, 0.0 } })

   tests.Shocktube3 = make_2state_setup
   ({ Pl = { 1.0, 1e+3, 0.0, 0.0, 0.0 },
      Pr = { 1.0, 1e-2, 0.0, 0.0, 0.0 } })

   tests.Shocktube4 = make_2state_setup
   ({ Pl = { 1.0, 1e-2, 0.0, 0.0, 0.0 },
      Pr = { 1.0, 1e+2, 0.0, 0.0, 0.0 } })

   tests.Shocktube5 = make_2state_setup
   ({ Pl = { 5.99924, 460.894, 19.59750, 0.0, 0.0 },
      Pr = { 5.99924,  46.095, -6.19633, 0.0, 0.0 } })

   tests.ContactWave = make_2state_setup
   ({ Pl = { 1.0, 1.0, 0.0, 0.7, 0.2 },
      Pr = { 0.1, 1.0, 0.0, 0.7, 0.2 } })

   tests.RMHDShocktube1 = make_2state_setup
   ({ Pl = { 1.000, 1.000, 0.000, 0.0, 0.0, 0.5, 1.0, 0.0 },
      Pr = { 0.125, 0.100, 0.000, 0.0, 0.0, 0.5,-1.0, 0.0 } })

   tests.RMHDShocktube2 = make_2state_setup
   ({ Pl = { 1.080, 0.950, 0.400, 0.3, 0.2, 2.0, 0.3, 0.3 },
      Pr = { 1.000, 1.000,-0.450,-0.2, 0.2, 2.5,-0.7, 0.5 } })

   tests.RMHDShocktube3 = make_2state_setup
   ({ Pl = { 1.000, 0.100, 0.999, 0.0, 0.0, 10.0, 0.7, 0.7 },
      Pr = { 1.000, 0.100,-0.999, 0.0, 0.0, 10.0,-0.7,-0.7 } })

   tests.RMHDShocktube4 = make_2state_setup
   ({ Pl = { 1.000, 5.000, 0.000, 0.3, 0.4, 1.0, 6.0, 2.0 },
      Pr = { 0.900, 5.300, 0.000, 0.0, 0.0, 1.0, 5.0, 2.0 } })

   tests.RMHDContactWave = make_2state_setup
   ({ Pl = { 10.0, 1.0, 0.0, 0.7, 0.2, 5.0, 1.0, 0.5 },
      Pr = {  1.0, 1.0, 0.0, 0.7, 0.2, 5.0, 1.0, 0.5 } })

   tests.RMHDRotationalWave = make_2state_setup
   ({ Pl = { 1, 1, 0.400000, -0.300000, 0.500000, 2.4, 1.00,-1.600000 },
      Pr = { 1, 1, 0.377347, -0.482389, 0.424190, 2.4,-0.10,-2.178213 } })

   return tests
