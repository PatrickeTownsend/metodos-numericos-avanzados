#-----Patrick Townsend-----#
using LinearAlgebra
using Plots
J(x) = [2*x[1] 2*x[2]; x[2] x[1]]

F(x) = [x[1]^2+x[2]^2-25; 
        x[1]*x[2]-9]

x = [10,0]
max_tol = 1e-8
max_iter = 5000
δ=1
iter = 0
Jac = J(x)
Fun = F(x)
while norm(δ)≥ max_tol
    iter +=1
    δ = -Jac\Fun
    x += δ

    println("iter: $iter  x: $x")
    Jac = J(x)
    Fun = F(x)
end
xf = x[1]
yf= x[2]
println("--------------------------------------------------")
println("converged at $iter iterations with x=$xf and y=$yf")
println("--------------------------------------------------")








