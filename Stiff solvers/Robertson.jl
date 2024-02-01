using Plots
using LinearAlgebra
using DifferentialEquations
using LaTeXStrings
include("Problems.jl")
include("Adams.jl")
include("RungeKutta.jl")

import .Problems
import .Adams
import .RungeKutta

IVP = Problems.IVP1
k = 3
h = 1e-6

@time t_1, y_1, f_1, N_1 = Adams.Bashforth(IVP, 4, h)
println("AB4 done! NFE: $N_1")
@time t_2, y_2, N_2 = RungeKutta.RK4(IVP,h)
println("RK4 done! NFE: $N_2")
@time t_3, y_3, N_3, J_3 = Adams.Moulton(IVP, k, 1e-4, 2)
println("AM3 done! NFE: $N_3 NJE: $J_3")
@time t_4, y_4, N_4, J_4 = RungeKutta.SDIRK4(IVP,3.4e-6)
println("SDIRK done! NFE:$N_4 NJE:$J_4")

#---------Exact solution aproximation---------#
f(u, p, t) = [-0.04*u[1] + (1e4)*u[2]*u[3]; 0.04*u[1] - (1e4)*u[2]*u[3] - (3e7)*u[2]^2; (3e7)*u[2]^2]
prob = ODEProblem(f, IVP.IC, IVP.tspan)
sol = solve(prob,KenCarp5())


#--------Error at t=tf------------#
error_1 = @. (norm(sol[end] - y_1[:,end])/norm(sol[end]))
error_2 = @. norm(sol[end] - y_2[:,end])/norm(sol[end])
error_3 = @. norm(sol[end] - y_3[:,end])/norm(sol[end])
error_4 = @. norm(sol[end] - y_4[:,end])/norm(sol[end])

#----------------Output-------------------#
println("------Results at t= tf-------------")
println("Adams-Bashforth 4 :", y_1[:,end], " error: ",norm(error_1))
println("Runge-Kutta 4 :", y_2[:,end], " error: ",norm(error_2))
println("Adams-Moulton 3 :", y_3[:,end], " error: ",norm(error_3))
println("SDIRK 4 :", y_4[:,end], " error: ",norm(error_4))
println("-----------------------------------")

#-----------Plotting--------------#
plot(t_3,y_3[3,:], title="Robertson problem", label="Adams - Bashforth 4",framestyle=:box, fontfamily="Computer Modern",
tickfontsize=11,legendfontsize=11,titlefontsize=14,labelfontsize=12,lw=2, legend=:bottomright)
plot!(t_3, y_3[3,:], label="Runge-Kutta 4")
plot!(t_3, y_3[3,:], label="Adams-Moulton 3")
p1 = plot!(t_4, y_4[3,:], label="SDIRK 4")
xlabel!("Time")
ylabel!(L"y_3(t)")
#savefig("Stiff solvers/plots/Robertson_3.pdf")