local tests = { }

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
               local r2 = x*x + y*y + z*z
               if r2 < 0.005 then
                  return { 1.000, 1.0, 0, 0, 0 }
               else
                  return { 0.125, 0.1, 0, 0, 0 }
               end
            end
         return g
      end
}


tests.ExplosionRmhd = {
   get_pinit =
      function(self, time)
         local g = 
            function(x,y,z)
               local r2 = x*x + y*y
               if r2 < 0.025 then
                  return { 1.000, 1.0, 0, 0, 0, 4, 0, 0 }
               else
                  return { 0.125, 0.1, 0, 0, 0, 4, 0, 0 }
               end
            end
         return g
      end
}

tests.KelvinHelmoltz = {
   get_pinit =
      function(self, time)
         local g = 
            function(x,y,z)
            local rho,vx,vy   
               if math.abs(y) > 0.25 then
                    rho = 1.0
                    vx = -0.5
                 else
                    rho = 2.0
                    vx =  0.5
                 end
            vx = 0.02*(math.random() - 0.5) + vx
            vy = 0.02*(math.random() - 0.5)
            return { rho, 2.5, vx, vy, 0.0 }   
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



tests.SrhdCase1_DFIM98 = {
   Pl = { 10.0, 13.30, 0.0, 0.0, 0.0 },
   Pr = {  1.0,  1e-6, 0.0, 0.0, 0.0 }
}

tests.SrhdCase2_DFIM98 = {
   Pl = { 1, 1e+3, 0.0, 0.0, 0.0 },
   Pr = { 1, 1e-2, 0.0, 0.0, 0.0 }
}

tests.MakeShocktubeProblem =
   function(states, opts)
      return {
	 get_pinit =
	    function(self)
	       local g =
		  function(x,y,z)
		     if not (opts and opts.reverse) then
			if x < 0.5 then
			   return states.Pl
			else
			   return states.Pr
			end
		     else
			if x < 0.5 then
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


return tests
