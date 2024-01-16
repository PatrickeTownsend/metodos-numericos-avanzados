using Plots
using LinearAlgebra

include("Problems.jl")
include("Adams.jl")
include("RungeKutta.jl")

import .Problems
import .Adams
import .RungeKutta

IVP = Problems.IVP2
k = 4
h = 1e-6

#t_1, y_1, f_1 = Adams.Bashforth(IVP, k, h)
t_2, y_2 = RungeKutta.RK4(IVP,h)
println(y_2[:,end])