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
h = 1e-6 #conervge a medias xd (3.4e-4)

t_1, y_1, f_1 = Adams.Bashforth(IVP, k, 1e-7)
println("AB4 done!")
t_2, y_2 = RungeKutta.RK4(IVP,h)
println("RK4 done!")
t_3, y_3 = Adams.Moulton(IVP, 3, h, 2) #0.7741664454953029 (1e-6)
println("AM3 done!")
t_4, y_4 = RungeKutta.SDIRK4(IVP,5e-6) # 0.7739429635005897 (5e-5)

#---------------#
println("------Results at t= tf-------------")
println("Adams-Bashforth 4 :", y_1[:,end])
println("Runge-Kutta 4 :", y_2[:,end])
println("Adams-Moulton 3 :", y_3[:,end])
println("SDIRK 4 :", y_4[:,end])
println("-----------------------------------")
plot(t_1,y_1', title="The Prothero -Robinson problem", label="Adams - Bashforth 4")
plot!(t_2, y_2', label="Runge-Kutta 4")
plot!(t_3, y_3', label="Adams-Moulton 3")
plot!(t_4, y_4', label="SDIRK 4")
xlabel!("Time")
ylabel!("y")