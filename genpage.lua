
local oo = require 'class'
local problems = require 'problems'

local vars = {'\\rho', 'p', 'v_x', 'v_y', 'v_z'}
for _,v in pairs(problems) do
   if oo.issubclass(v, problems.TwoStateProblem) then
      print('* '..oo.classname(v))
      print('\n|---+---+---|')
      for i=1,5 do
	 print(string.format('| $%s$ | %f | %f|', vars[i], v.state1[i], v.state2[i]))
      end
      print('|---+---+---|\n')
   end
end
