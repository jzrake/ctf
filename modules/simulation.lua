
local oo = require 'class'

local SimulationBase = oo.class('SimulationBase')

function SimulationBase:_notimplemented()
   error("method "..oo.classname(self)..':'..debug.getinfo(2, 'n').name..
	 " needs to be implemented")
end
function SimulationBase:initialize_behavior() end
function SimulationBase:initialize_solver() self:_notimplemented() end
function SimulationBase:initialize_physics() self:_notimplemented() end
function SimulationBase:finalize_solver() self:_notimplemented() end
function SimulationBase:advance_physics() self:_notimplemented() end
function SimulationBase:set_time_increment() self:_notimplemented() end
function SimulationBase:checkpoint_condition()
   return (self.status.simulation_time - self.status.last_checkpoint_time >=
	   self.behavior.checkpoint_cadence)
end
function SimulationBase:checkpoint_write() end
function SimulationBase:checkpoint_read() end
function SimulationBase:user_work_iteration() end
function SimulationBase:user_work_finish() end
function SimulationBase:advance_status()
   local status = self.status
   status.iteration_number = status.iteration_number + 1
   status.simulation_time = status.simulation_time + status.time_increment
end
function SimulationBase:iteration_message()
   local n = self.status.iteration_number
   if n % self.behavior.message_cadence ~= 0 then return end
   local msg = string.format("%05d: t=%3.2f dt=%2.1e",
			     self.status.iteration_number,
			     self.status.simulation_time,
			     self.status.time_increment)
   print(msg)
end
function SimulationBase:continue_condition()
   return
      self.status.iteration_number < self.behavior.max_iteration and
      self.status.wall_runtime < self.behavior.max_wall_runtime and
      self.status.simulation_time < self.behavior.max_simulation_time
end
function SimulationBase:main_loop()
   while self:continue_condition() do
      if self:checkpoint_condition() then
	 self:checkpoint_write()
	 self.status.checkpoint_number = self.status.checkpoint_number + 1
	 self.status.last_checkpoint_time = self.status.simulation_time
      end
      self:user_work_iteration()
      self:set_time_increment()
      self:advance_physics()
      self:advance_status()
      self:iteration_message()
   end
end
function SimulationBase:run(problem)
   self.problem = problem
   problem.simulation = self
   self:initialize_behavior()
   self:initialize_solver()
   self:initialize_physics()
   self:main_loop()
   self:user_work_finish()
   self:finalize_solver()
end
function SimulationBase:__init__(user_opts)
   self.status   = { iteration_number     = 0,
		     wall_runtime         = 0.0,
		     simulation_time      = 0.0,
		     time_incremement     = 0.0,
		     checkpoint_number    = 0,
		     last_checkpoint_time = 0.0 }
   self.behavior = { max_iteration        = math.huge,
		     max_wall_runtime     = math.huge,
		     max_simulation_time  = 1.0, -- in simulation time
		     checkpoint_cadence   = 1.0, -- in simulation time
		     message_cadence      = 1 }  -- in iterations
   self.user_opts = user_opts
end


local function main()
   local optparse = require 'optparse'
   local parser = optparse.OptionParser{usage="%prog [options] [input_args]",
					version="CTF version 1.0"}
   parser.add_option{"-s", "--solver", dest="solver", help="solver type"}

   local opts, args = parser.parse_args()

   local MySimulation = oo.class('MySimulation', SimulationBase)
   function MySimulation:initialize_physics()
      print'init physics'
   end
   function MySimulation:set_time_increment()
      self.status.time_increment = 0.01
   end
   function MySimulation:advance_physics()
      
   end
   local sim = MySimulation()
   
   sim.initialize_solver = function() print 'init solver' end

   sim:initialize_behavior()
   sim:initialize_solver()
   sim:initialize_physics()
   sim:main_loop()
end

if ... then -- if __name__ == "__main__"
   return {SimulationBase=SimulationBase}
else
   main()
   print(debug.getinfo(1).source, ": All tests passed")
end
