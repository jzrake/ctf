
local array  = require 'array'
local visual = require 'visual'

local Nx = 400
local Ny = 200
local A = array.array{Nx, Ny, 2}
local Avec = A:vector()
for i=0,Nx-1 do
   for j=0,Ny-1 do
      Avec[(j * Nx + i) * 2 + 0] = -(j-200)/Ny
      Avec[(j * Nx + i) * 2 + 1] =  (i-200)/Nx
   end
end
visual.line_integral_convolution(A:buffer(), Nx, Ny, "lic.ppm")
