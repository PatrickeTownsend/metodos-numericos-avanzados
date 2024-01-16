using Plots
using LinearAlgebra

include("Problems.jl")
include("Adams.jl")
include("RungeKutta.jl")

import .Problems
import .Adams
import .RungeKutta

IVP = Problems.IVP1
k = 4
h = 1e-3

#t_1, y_1, f_1 = Adams.Bashforth(IVP, k, h)
#t_2, y_2 = RungeKutta.RK4(IVP,h)
t_3, y_3 = Adams.Moulton(IVP, k, h, 1)
#println(y_3[:,end])
#println(y_2[:,end])
plot(t_3,y_3[1,:])