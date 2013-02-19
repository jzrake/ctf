
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
local MyMara = oo.class('MyMara', MyBase, MaraSim.MaraSimulation)
local MyFish = oo.class('MyFish', MyBase, FishSim.FishSimulation)


function MyFish:initialize_solver()
   oo.super(self, FishSim.FishSimulation):initialize_solver()
   self.Primitive = self.grid.primitive
   self.Gravity = self.grid.gravity
   self.dx = 1.0 / self.N
end

function MyBase:_map_solution(t)
   local N  = self.N
   local Ng = self.Ng
   local dx = 1.0 / N

   local P = array.array{N + 2*Ng, 5}
   local Pvec = P:vector()

   for n=0,#Pvec/5-1 do
      local x  = (n - Ng + 0.5) * dx
      local Pi = self.problem:solution(x,0,0,t)
      Pvec[5*n + 0] = Pi[1]
      Pvec[5*n + 1] = Pi[2]
      Pvec[5*n + 2] = Pi[3]
      Pvec[5*n + 3] = Pi[4]
      Pvec[5*n + 4] = Pi[5]
   end
   return P
end

function MyBase:initialize_behavior()
   local opts = self.user_opts
   local cpi = opts.cpi or 1.0
   local tmax = opts.tmax or self.problem:finish_time()
   self.behavior.message_cadence = opts.message_cadence or 10
   self.behavior.checkpoint_cadence = tonumber(cpi)
   self.behavior.max_simulation_time = tonumber(tmax)
end

function MyBase:checkpoint_write(fname)
   if not next(hdf5) then -- next(t) is true when t is empty
      print('warning! Could not write checkpoint, HDF5 not available')
      return
   end

   local base = self.user_opts.id or 'chkpt'
   local n = self.status.checkpoint_number
   local t = self.status.simulation_time
   local Pexact = self:_map_solution(t)
   local Ng = self.Ng
   local fname = fname or string.format('data/%s.%04d.h5', base, n)

   os.execute('mkdir -p data')
   print('writing checkpoint ' .. fname)

   local outfile = hdf5.File(fname, 'w')
   outfile['prim' ]   = self.Primitive[{{Ng,-Ng},nil}]
   outfile['exact']   = Pexact        [{{Ng,-Ng},nil}]
   outfile['id']      = base
   outfile['time']    = t
   outfile['problem'] = oo.classname(self.problem)
   outfile:close()
end

function MyBase:user_work_finish()
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

   if self.user_opts.plot and not self.user_opts.convergence then
      util.plot{['code' ]=P     [{{Ng,-Ng},{0,1}}]:table(),
                ['exact']=Pexact[{{Ng,-Ng},{0,1}}]:table()}
   end

   local output = self.user_opts.output

   if output and not self.user_opts.convergence then
      if util.endswith(output, '.h5') then
         self:checkpoint_write(output)
      else
         local f = io.open(output, 'w')
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

   self.problem:user_work_finish()
end

function MyBase:user_work_iteration()
   self.problem:user_work_iteration()
end



local function main()
   local usage = "tests-1d <problem> [<options>]"
   local parser = optparse.OptionParser{usage=usage,
                                        version="CTF version 1.0"}

   parser.add_option{"--cpi", dest="cpi", help="checkpoint interval"}
   parser.add_option{"--id", dest="id", help="problem ID: used for checkpoint names"}
   parser.add_option{"--cfl", dest="CFL",
                     help="Courant-Freidrichs-Lewy time-step constraint"}
   parser.add_option{"--tmax", dest="tmax", help="end simulation time"}
   parser.add_option{"--plot", dest="plot", action="store_true"}
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

      --util.plot({['L1 error']=ErrorTable}, {cmds={'set logscale'}})

      if opts.output then
         local f = io.open(opts.output, 'w')
         for N,L in pairs(ErrorTable) do
            f:write(string.format('%d %8.6e\n', N, L))
         end
         f:close()
      end
   else
      local sim = sim_class(opts)
      local problem = problem_class(opts)
      sim:run(problem)
   end
end

main()
