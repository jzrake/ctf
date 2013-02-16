
local FishCls = require 'FishClasses'
local D = FishCls.FluidDescriptor{gamma=1.4, fluid='nrhyd'}
local R = FishCls.RiemannSolver()
local SL = FishCls.FluidState(D)
local SR = FishCls.FluidState(D)

SL.primitive[0] = 1.000 -- density
SR.primitive[0] = 0.125
SL.primitive[1] = 1.0   -- pressure
SR.primitive[1] = 0.1
SL.primitive[2] = 0.0   -- velocity
SR.primitive[2] = 0.0

local x = 1.0
local t = 0.5
local P = R:solve(SL, SR, x/t)
print(P)
