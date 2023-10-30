using Plots
using LinearAlgebra
using Random


include("Functions.jl")
N = 32 # Number of charges
max_iter = 10000 # Max iteration 
max_tol = 1e-10 # Convergence tolerance
N_iter = 1
r = zeros(N,3)
U = zeros(max_iter)
α = 0.01
r[:,:] = Initialization(N)
U[1] = PotentialEnergy(r,N)
#r[:,:] = InitRandom(N)
ϵ= zeros(N,3,max_iter)
while N_iter < max_iter 
    tot_grad = zeros(3)
    for i = 1:N
        r[i, :] = r[i,:] - α*Gradient(r,i,N)
        r[i, :] /= norm(r[i,:])
        tot_grad += Gradient(r,i,N)
    end
    U[N_iter+1] += PotentialEnergy(r,N)
    potential = 
    println("Iteration $N_iter Total Grad ",norm(tot_grad)," x: ",abs(tot_grad[1])," y: ",abs(tot_grad[2]), " z: ", abs(tot_grad[3]))
    #println("Iter $N_iter - ",PotentialEnergy(r,N))
    # if abs(U[N_iter+1]-U[N_iter]) ≤ max_tol
    #      break
    # end
    if norm(tot_grad) ≤ max_tol
        break
   end

    N_iter += 1
end

U = PotentialEnergy(r,N)


x = r[:,1]
y = r[:,2]
z = r[:,3]
PlotSphere(x,y,z,"Gradient")
println("-----Gradient Descent-----")
println("Potential: $U")
println(r[1,:])
println(r[2,:])



        