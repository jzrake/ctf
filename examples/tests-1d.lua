
local optparse = require 'optparse'
local FishSim  = require 'FishSim'
local MaraSim  = require 'MaraSim'
local problems = require 'problems'
local util     = require 'util'

local cmds = { }

-- function cmds.convergence(sim_class, problem, opts)
--    local ErrorTable = { }
--    for _,res in pairs{8,16,32,64,128,256,512,1024} do
--       opts.N = res
--       local sim = sim_class(opts)
--       sim:run(problem)
--       ErrorTable[sim.N] = sim.L1error
--    end
--    util.pretty_print(ErrorTable)
--    print("estimated convergence rate is",
-- 	 math.log10(ErrorTable[512]/ErrorTable[64]) / math.log10(512/64))
--    util.plot({['L1 error']=ErrorTable}, {'set logscale'})
-- end

local function main()
   local usage = "test-1d <problem> [<options>]"
   local parser = optparse.OptionParser{usage=usage,
					version="CTF version 1.0"}

   parser.add_option{"--cpi", dest="cpi", help="checkpoint interval"}
   parser.add_option{"--cfl", dest="CFL",
		     help="Courant-Freidrichs-Lewy time-step constrain"}
   parser.add_option{"--tmax", dest="tmax", help="end simulation time"}
   parser.add_option{"--plot", dest="plot", action="store_true"}
   parser.add_option{"--problem", dest="problem", help="problem name to run"}
   parser.add_option{"--reconstruction", dest="reconstruction"}
   parser.add_option{"--riemann", dest="riemann"}
   parser.add_option{"--advance", dest="advance",
		     help="which Runge-Kutta to use for solution advance"}
   parser.add_option{"--solver", dest="solver"}
   parser.add_option{"--backend", "-b", dest="backend",
		     help="which backend code to use: currently Mara or Fish"}
   parser.add_option{"--resolution, -N", dest="resolution", help="grid resolution"}
   parser.add_option{"--self-gravity", dest="self_gravity", action="store_true",
		     help="include self gravity"}
   parser.add_option{"--convergence", dest="convergence", action="store_true",
		     help="run a convergence test over a few resolutions"}

   local opts, args = parser.parse_args()
   local problem_class = problems[arg[2]]

   if not problem_class then
      print('usage: '..usage)
      print("valid problem names are:")
      util.pretty_print(problems, '\t')
      return
   end

   local sim_class = (
      {mara=MaraSim.MaraSimulation,
       fish=FishSim.FishSimulation})[(opts.backend or 'mara'):lower()]

   if not sim_class then
      print("valid backend codes are:")
      util.pretty_print({'Mara', 'Fish'}, '\t')
      return
   end

   local sim = sim_class(opts)
   local problem = problem_class(opts)
   sim:run(problem)
end

main()
