using Plots
using LinearAlgebra
using Random
include("Functions.jl")

#---Set up---#
N = 32 # Number of charges
max_iter = 10000 # Max iteration 
max_tol = 1e-18 # Convergence tolerance
N_iter = 1 #Iteration counter
α = 0.01 #Learning rate

#---Arrays---#
r = zeros(N,3)
U = zeros(2)

#---Initial Guess---#
#r[:,:] = Initialization(N)
r[:,:] = InitRandom(N)
U[1] = PotentialEnergy(r,N)

@time while N_iter < max_iter 
    tot_grad = zeros(3)
    for i = 1:N
        r[i, :] = r[i,:] - α*Gradient(r,i,N)
        r[i, :] /= norm(r[i,:])
        tot_grad += Gradient(r,i,N)
    end
    U[2] = PotentialEnergy(r,N)
    println("Iter $N_iter) U: ",PotentialEnergy(r,N)," x: ",abs(tot_grad[1])," y: ",abs(tot_grad[2]), " z: ", abs(tot_grad[3]))
    if abs(U[2]-U[1]) ≤ max_tol
        break
    end
    N_iter += 1
    U[1] = U[2]
end

#---Particles position---#
x = r[:,1]
y = r[:,2]
z = r[:,3]

#---Output---#
println("-----Gradient Descent-----")
println("Converged at iter $N_iter with U = ",U[2])

#---Plot---#
PlotSphere(x,y,z,"Gradient")





        