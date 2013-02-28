
local array  = require 'array'
local visual = require 'visual'

local Nx = 400
local Ny = 200
local A = array.array{Nx, Ny, 2}
local B = array.array{Nx, Ny}
local Avec = A:vector()
local Bvec = B:vector()
for i=0,Nx-1 do
   for j=0,Ny-1 do
      local x = (i-200)/Nx
      local y = (j-100)/Ny
      local m = i*Ny + j
      Avec[2*m + 0] = -y
      Avec[2*m + 1] =  x
      Bvec[m] = math.sin(4 * math.pi * x)
   end
end
visual.line_integral_convolution(A:buffer(), Nx, Ny, "lic.ppm")
visual.write_ppm(B:buffer(), Nx, Ny, "img.ppm")

os.execute('convert lic.ppm lic.png')
os.execute('convert img.ppm img.png')
