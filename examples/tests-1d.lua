
local optparse = require 'optparse'
local FishSim  = require 'FishSim'
local problems = require 'problems'
local util     = require 'util'

local cmds = { }

function cmds.runandplot(problem, opts)
   local sim = FishSim.FishSimulation(opts)
   sim:run(problem)
end

function cmds.convergence(problem, opts)
   local ErrorTable = { }
   for _,res in pairs{8,16,32,64,128,256,512,1024} do
      opts.N = res
      local sim = FishSim.FishSimulation(opts)
      sim:run(problem)
      ErrorTable[sim.N] = sim.L1error
   end
   util.pretty_print(ErrorTable)
   print("estimated convergence rate is",
	 math.log10(ErrorTable[512]/ErrorTable[64]) / math.log10(512/64))
   util.plot({['L1 error']=ErrorTable}, {'set logscale'})
end

local function main()
   local parser = optparse.OptionParser{usage="%prog [options] [input_args]",
					version="CTF version 1.0"}

   parser.add_option{"--cpi", dest="cpi", help="checkpoint interval"}
   parser.add_option{"--tmax", dest="tmax", help="end simulation time"}
   parser.add_option{"--plot", dest="plot", action="store_true"}
   parser.add_option{"--problem", dest="problem", help="problem name to run"}
   parser.add_option{"--reconstruction", dest="reconstruction"}
   parser.add_option{"--riemann", dest="riemann"}
   parser.add_option{"--advance", dest="advance"}
   parser.add_option{"--solver", dest="solver"}
   parser.add_option{"-N", dest="N", help="resolution"}

   local opts, args = parser.parse_args()

   local problem_class = problems[opts.problem or 'densitywave']
   if not problem_class then
      print("valid problem names are:")
      util.pretty_print(problems, '\t')
      return
   end
   if not cmds[args[2]] then
      print("valid sub commands are:")
      util.pretty_print(cmds, '\t')
   else
      cmds[args[2]](problem_class(opts), opts)
   end
end

main()
