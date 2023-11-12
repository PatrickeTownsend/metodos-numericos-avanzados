using LinearAlgebra
include("Functions.jl")
function gradient_descent(N::Int, r₀::Array, max_tol::Float64, max_iter::Int,α::Float64, type::String)
    r = zeros(N,3)
    U = zeros(2)
    r[:,:] = r₀
    U[1] = PotentialEnergy(r,N)
    N_iter = 1
    iterations = zeros(Int64,max_iter)
    derivatives = zeros(max_iter,3)
    residuals = zeros(max_iter)

    while N_iter < max_iter 
        tot_grad = zeros(3)
        for i = 1:N
            r[i, :] = r[i,:] - α*Gradient(r,i,N)
            r[i, :] /= norm(r[i,:])
            tot_grad += Gradient(r,i,N)
        end
        U[2] = PotentialEnergy(r,N)
        derivatives[N_iter,:] += abs.(tot_grad)
        residuals[N_iter] += abs((U[2]-U[1]))
        if abs(U[2]-U[1]) ≤ max_tol
            break
        end
        N_iter += 1
        iterations[N_iter]+=N_iter
        U[1] = U[2]
    end
    derivatives = derivatives[1:N_iter,:]
    iterations = iterations[1:N_iter]
    residuals = residuals[1:N_iter]
    grad = norm(derivatives[N_iter,:])
    #---Output---#
    println("                                  ")
    println("----------------------------------")
    println("---------Gradient Descent---------")
    println("----------------------------------")
    println("Initial estimation: $type")
    println("Converged at iter $N_iter with U = ",U[2]," Residual: ",grad)
    println("                                        ")
    println("Iter $N_iter) x: ",derivatives[N_iter,1]," y: ",derivatives[N_iter,2]," z: ",derivatives[N_iter,3])
    println("--------------------Execution Time----------------------------")
    return r,residuals,iterations,derivatives
end





        