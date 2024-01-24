using Plots
using LinearAlgebra

include("Problems.jl")
include("Adams.jl")
include("RungeKutta.jl")

import .Problems
import .Adams
import .RungeKutta

IVP = Problems.IVP1
k = 3
h = 1e-6 #conervge a medias xd (3.4e-4)

t_1, y_1, f_1 = Adams.Bashforth(IVP, 4, h)
println("AB4 done!")
t_2, y_2 = RungeKutta.RK4(IVP,h)
println("RK4 done!")
t_3, y_3 = Adams.Moulton(IVP, k, 1e-5, 2)
println("AM3 done!")
t_4, y_4 = RungeKutta.SDIRK4(IVP,3.4e-4)

#---------------#
println("------Results at t= tf-------------")
println("Adams-Bashforth 4 :", y_1[:,end])
println("Runge-Kutta 4 :", y_2[:,end])
println("Adams-Moulton 3 :", y_3[:,end])
println("SDIRK 4 :", y_4[:,end])
println("-----------------------------------")

plot(t_1,y_1[1,:], title="Robertson problem", label="Adams - Bashforth 4")
plot!(t_2, y_2[1,:], label="Runge-Kutta 4")
plot!(t_3, y_3[1,:], label="Adams-Moulton 3")
plot!(t_4, y_4[1,:], label="SDIRK 4")
xlabel!("Time")
ylabel!("y₁")

