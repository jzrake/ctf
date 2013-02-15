
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
function SimulationBase:local_mesh_size() self:_notimplemented() end
function SimulationBase:set_time_increment() self:_notimplemented() end
function SimulationBase:checkpoint_condition()
   return (self.status.simulation_time - self.status.last_checkpoint_time >=
	   self.behavior.checkpoint_cadence)
end
function SimulationBase:checkpoint_write() end
function SimulationBase:checkpoint_read() end
function SimulationBase:user_work_iteration() end
function SimulationBase:user_work_finish() end
function SimulationBase:report_configuration() end
function SimulationBase:advance_status()
   local status = self.status
   local kz = self:local_mesh_size() / 1000
   status.iteration_number = status.iteration_number + 1
   status.simulation_time = status.simulation_time + status.time_increment
   status.kilozones_per_second = kz / self.profiler['advance_physics']
end
function SimulationBase:iteration_message()
   local n = self.status.iteration_number
   if n % self.behavior.message_cadence ~= 0 then return end
   local msg = string.format("%05d: t=%3.2f dt=%2.1e %3.2fkz/s",
			     self.status.iteration_number,
			     self.status.simulation_time,
			     self.status.time_increment,
			     self.status.kilozones_per_second)
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
      self:time('user_work_iteration')
      self:time('set_time_increment')
      self:time('advance_physics')
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
   self:report_configuration()
   self:main_loop()
   self:user_work_finish()
   self:finalize_solver()
   print('\nProfiler output:')
   print(string.format("%30s     %s", 'function name', 'cumulative time'))
   print(string.format("%30s     %s", '-------------', '---------------'))
   for k,v in pairs(self.profiler) do
      print(string.format("%30s ... %3.2es", k,v))
   end
end
function SimulationBase:__init__(user_opts)
   self.status   = { iteration_number     = 0,
		     wall_runtime         = 0.0,
		     simulation_time      = 0.0,
		     time_incremement     = 0.0,
		     checkpoint_number    = 0,
		     last_checkpoint_time = 0.0,
		     kilozones_per_second = 0.0 }
   self.behavior = { max_iteration        = math.huge,
		     max_wall_runtime     = math.huge,
		     max_simulation_time  = 1.0, -- in simulation time
		     checkpoint_cadence   = 1.0, -- in simulation time
		     message_cadence      = 1 }  -- in iterations
   self.profiler = { }
   self.user_opts = user_opts or { }
end
function SimulationBase:time(func_name, ...)
   local start = os.clock()
   self[func_name](self, ...)
   self.profiler[func_name] = os.clock() - start
end

return {SimulationBase=SimulationBase}
