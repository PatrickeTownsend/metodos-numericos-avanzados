using Plots
using LinearAlgebra
include("Functions.jl")
N = 32 # Number of charges
max_iter = 10000 # Max iteration 
max_tol = 1e-8 # Convergence tolerance
N_iter = 1
r = zeros(N,3)
α = 0.01

r[:,:] = Initialization(N)
while N_iter < max_iter 
    tot_grad = zeros(3)
    for i = 1:N
        r[i, :] = r[i,:] - α*Gradient(r,i,N)
        r[i, :] /= norm(r[i,:])
        tot_grad += Gradient(r,i,N)
    end

    if norm(tot_grad) ≤ max_tol
        break
    end

    N_iter += 1

end

U = PotentialEnergy(r,N)


x = r[:,1]
y = r[:,2]
z = r[:,3]
#PlotSphere(x,y,z)
println("-----Gradient Descent-----")
println("Potential: $U")
println(r[1,:])
println(r[2,:])
        


        