using LinearAlgebra
include("Functions.jl")
function gradient_descent(N::Int, r₀::Array, max_tol::Float64, max_iter::Int,α::Float64)
    r = zeros(N,3)
    U = zeros(2)
    r[:,:] = r₀
    U[1] = PotentialEnergy(r,N)
    N_iter = 1
    iterations = zeros(Int64,max_iter)
    residuals = zeros(max_iter,3)
    while N_iter < max_iter 
        tot_grad = zeros(3)
        for i = 1:N
           r[i, :] = r[i,:] - α*Gradient(r,i,N)
           r[i, :] /= norm(r[i,:])
           tot_grad += Gradient(r,i,N)
        end
        U[2] = PotentialEnergy(r,N)
        residuals[N_iter,:] += abs.(tot_grad)
        if abs(U[2]-U[1]) ≤ max_tol
            break
        end
        N_iter += 1
        iterations[N_iter]+=N_iter
        U[1] = U[2]
    end
    residuals = residuals[1:N_iter,:]
    iterations = iterations[1:N_iter]
    #---Output---#
    println("--------------------------")
    println("-----Gradient Descent-----")
    println("--------------------------")
    println("Converged at iter $N_iter with U = ",U[2])
    PlotSphere(r[:,1],r[:,2],r[:,3],"Gradient")
    PlotResiduals(residuals,iterations,"Gradient")
    return nothing
end





        