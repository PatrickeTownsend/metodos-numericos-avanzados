using Plots
using LinearAlgebra
using DifferentialEquations
include("Problems.jl")
include("Adams.jl")
include("RungeKutta.jl")

import .Problems
import .Adams
import .RungeKutta

IVP = Problems.IVP2
k = 4


@time t_1, y_1, f_1, N_1 = Adams.Bashforth(IVP, k, 1e-7)
println("AB4 done! NFE: $N_1")
@time t_2, y_2, N_2 = RungeKutta.RK4(IVP,1e-7)
println("RK4 done! NFE: $N_2")
@time t_3, y_3, N_3, J_3 = Adams.Moulton(IVP, 3, 1e-7, 2)
println("AM3 done! NFE: $N_3 NJE: $J_3")
@time t_4, y_4, N_4, J_4 = RungeKutta.SDIRK4(IVP,2.35e-6)
println("SDIRK done! NFE:$N_4 NJE:$J_4")

#------Error with AB5--------#
t_1_e, y_1_e, f_1_e, N_1_e = Adams.Bashforth(IVP,5,1e-7)
t_2_e, y_2_e, f_2_e, N_2_e = Adams.Bashforth(IVP,5,1e-7)
t_4_e, y_4_e, f_3_e, N_4_e = Adams.Bashforth(IVP,5,2.35e-7)

e1 = @. abs(y_1_e - y_1)/y_1_e
e2 = @. abs(y_2_e - y_2)/y_2_e
e3 = @. abs(y_2_e - y_3)/y_2_e
e4 = @. abs(y_4_e - y_4)/y_4_e

#------Exact solution aproximation------------#
f(u, p, t) = @. cos(π/4 + t) - (1e6)*(u - sin(π/4 + t))
prob = ODEProblem(f, IVP.IC, IVP.tspan)
sol = solve(prob,KenCarp5())

#--------Error at t=tf------------#
error_1 = @. (abs(sol[4] - y_1[end])/abs(sol[4]))
error_2 = @. abs(sol[4] - y_2[end])/abs(sol[4])
error_3 = @. abs(sol[4] - y_3[end])/abs(sol[4])
error_4 = @. abs(sol[4] - y_4[end])/abs(sol[4])

#----------------------Output--------------------------------#
println("------Results at t= tf-------------")
println("Adams-Bashforth 4 :", y_1[:,end], " error: ",error_1)
println("Runge-Kutta 4 :", y_2[:,end], " error: ",error_2)
println("Adams-Moulton 3 :", y_3[:,end], " error: ",error_3)
println("SDIRK 4 :", y_4[:,end], " error: ",error_4)
println("-----------------------------------")

plot(t_1,(e1'), title="The Prothero - Robinson problem", label="Adams - Bashforth 4",framestyle=:box, fontfamily="Computer Modern",
tickfontsize=11,legendfontsize=11,titlefontsize=14,labelfontsize=12,lw=2, grid=true)
plot!(t_2, (e2'), label="Runge-Kutta 4",lw=2)
plot!(t_3, (e3'), label="Adams-Moulton 3",lw=2, legend=:topright)
plot!(t_4, e4', label="SDIRK 4",lw=2)
xlabel!("Time")
p1 = ylabel!("error relativo")
#savefig("Stiff solvers/plots/protero_error.pdf")
display(p1)

#------------Plotting------------#
plot(t_1,y_1', title="The Prothero - Robinson problem", label="Adams - Bashforth 4",framestyle=:box, fontfamily="Computer Modern",
tickfontsize=11,legendfontsize=11,titlefontsize=14,labelfontsize=12,lw=2, grid=true)
plot!(t_2, y_2', label="Runge-Kutta 4",lw=2)
plot!(t_3, y_3', label="Adams-Moulton 3",lw=2)
p = plot!(t_4, y_4', label="SDIRK 4",lw=2)
xlabel!("Time")
ylabel!("y(t)")
#savefig("Stiff solvers/plots/protero.pdf")
display(p)