
local optparse = require 'optparse'
local FishSim  = require 'FishSim'
local MaraSim  = require 'MaraSim'
local problems = require 'problems'
local util     = require 'util'
local oo       = require 'class'
local sim      = require 'simulation'
local array    = require 'array'
local hdf5     = require 'lua-hdf5.LuaHDF5'

local MyBase = oo.class('MyBase', sim.SimulationBase)
local MyMara = oo.class('MyMara', MaraSim.MaraSimulation)
local MyFish = oo.class('MyFish', FishSim.FishSimulation)

local function build_mysim(cls)

   function cls:_map_solution(t)
      local N  = self.N
      local Ng = self.Ng
      local dx = self.dx

      local P = array.array{N + 2*Ng, 5}
      local Pvec = P:vector()

      for n=0,#Pvec/5-1 do
         local x  = (n - Ng) * dx
         local Pi = self.problem:solution(x,0,0,t)
         Pvec[5*n + 0] = Pi[1]
         Pvec[5*n + 1] = Pi[2]
         Pvec[5*n + 2] = Pi[3]
         Pvec[5*n + 3] = Pi[4]
         Pvec[5*n + 4] = Pi[5]
      end
      return P
   end

   function cls:initialize_behavior()
      local opts = self.user_opts
      local cpi = opts.cpi or 1.0
      local tmax = opts.tmax or self.problem:finish_time()
      local dynamical_time = self.problem:dynamical_time()
      self.behavior.message_cadence = opts.message_cadence or 10
      self.behavior.checkpoint_cadence = cpi * dynamical_time
      self.behavior.max_simulation_time = tmax * dynamical_time
   end

   function cls:checkpoint_write()
      if not hdf5.File then
	 print('warning! Could not write checkpoint, HDF5 not available')
	 return
      end

      local n = self.status.checkpoint_number
      local t = self.status.simulation_time
      local Ng = self.Ng

      local Pexact = self:_map_solution(t)

      os.execute('mkdir -p data')

      local fname = string.format('data/chkpt.%04d.h5', n)
      local outfile = hdf5.File(fname, 'w')
      outfile['prim' ] = self.Primitive[{{Ng,-Ng},nil}]
      outfile['grav' ] = self.Gravity  [{{Ng,-Ng},nil}]
      outfile['exact'] = Pexact        [{{Ng,-Ng},nil}]
      outfile:close()
   end

   function cls:user_work_finish()
      local t  = self.status.simulation_time
      local P  = self.Primitive
      local Ng = self.Ng

      local Pexact = self:_map_solution(t)
      local dx = self.dx
      local P0 = Pexact[{{Ng,-Ng},nil}]:vector()
      local P1 = P     [{{Ng,-Ng},nil}]:vector()
      local L1 = 0.0

      for i=0,#P0/5-1,5 do
         L1 = L1 + math.abs(P1[i] - P0[i]) * dx
      end
      self.L1error = L1

      if self.user_opts.plot then
         util.plot{['code' ]=P     [{{Ng,-Ng},{0,1}}]:table(),
                   ['exact']=Pexact[{{Ng,-Ng},{0,1}}]:table()}
      end

      if self.user_opts.output then
         local f = io.open(self.user_opts.output, 'w')
         local P0 = P[{{Ng,-Ng},{0,1}}]:table()
         local P1 = P[{{Ng,-Ng},{1,2}}]:table()
         local P2 = P[{{Ng,-Ng},{2,3}}]:table()
         local P3 = P[{{Ng,-Ng},{3,4}}]:table()
         local P4 = P[{{Ng,-Ng},{4,5}}]:table()
         for i=1, self.N do
            local line = string.format(
               "%d %+12.8e %+12.8e %+12.8e %+12.8e %+12.8e\n",
               i, P0[i], P1[i], P2[i], P3[i], P4[i])
            f:write(line)
         end
         f:close()
      end
   end
end

build_mysim(MyMara)
build_mysim(MyFish)


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
   parser.add_option{"--resolution", "-N", dest="resolution", help="grid resolution"}
   parser.add_option{"--message-cadence", dest="message_cadence",
                     help="print a message every N iterations"}
   parser.add_option{"--self-gravity", dest="self_gravity", action="store_true",
                     help="include self gravity"}
   parser.add_option{"--convergence", dest="convergence", action="store_true",
                     help="run a convergence test over a few resolutions"}
   parser.add_option{"--output", "-o", dest="output",
                     help="write an ASCII table of the final solution"}

   local opts, args = parser.parse_args()
   local problem_class = problems[arg[2]]

   if not problem_class then
      print('usage: '..usage)
      print("valid problem names are:")
      util.pretty_print(problems.problems_1d, '\t')
      return
   end

   local sim_class = (
      {mara=MyMara,
       fish=MyFish})[(opts.backend or 'fish'):lower()]

   if not sim_class then
      print("valid backend codes are:")
      util.pretty_print({'Mara', 'Fish'}, '\t')
      return
   end

   if opts.convergence then
      local ErrorTable = { }
      for _,res in pairs{8,16,32,64,128,256,512,1024} do
         opts.resolution = res
         local sim = sim_class(opts)
         local problem = problem_class(opts)
         sim:run(problem)
         ErrorTable[res] = sim.L1error
      end
      util.pretty_print(ErrorTable)
      print("estimated convergence rate is",
	    math.log10(ErrorTable[128]/ErrorTable[64]) / math.log10(128/64))
      util.plot({['L1 error']=ErrorTable}, {cmds={'set logscale'}})
   else
      local sim = sim_class(opts)
      local problem = problem_class(opts)
      sim:run(problem)
   end
end

main()
