using LinearAlgebra
include("Functions.jl")
function newton2nd(N::Int, r₀::Array, max_tol::Float64, max_iter::Int,α::Float64)
    r = zeros(N,3,2)
    U = zeros(2)
    r[:,:,1] += r₀
    s = zeros(2,3)
    U[1] = PotentialEnergy(r[:,:,1],N)
    N_iter = 1
    iterations = zeros(Int64,max_iter)
    residuals = zeros(max_iter,3)
    while N_iter < max_iter
        tot_grad = zeros(3)
        for i = 1:N
            if N_iter == 1
                s[1,:] =  -α*Gradient(r[:,:,1],i,N)
                r[i,:,2] = r[i,:,1] + s[1,:]
                r[i,:,2] /= norm(r[i,:,2])
            else
                r[i,:,2] = r[i,:,1] + α*s[1,:]
                r[i,:,2] /= norm(r[i,:,2])
                s[2,:] = -Gradient(r[:,:,2],i,N) + s[1,:]*(Gradient(r[:,:,2],i,N)'*Gradient(r[:,:,2],i,N))/(Gradient(r[:,:,1],i,N)'*Gradient(r[:,:,1],i,N))
                tot_grad+=Gradient(r[:,:,2],i,N)
            end
        end
            U[2] = PotentialEnergy(r[:,:,2],N)
            residuals[N_iter,:]+=tot_grad
            if abs(U[2]-U[1]) ≤ max_tol
                break
            else
                N_iter+=1
                U[1]=U[2]
                r[:,:,1] = r[:,:,2]
                s[1,:] = s[2,:]
                iterations[N_iter]+=N_iter
            end
    end
    residuals = residuals[1:N_iter,:]
    iterations = iterations[1:N_iter]     
    #---Output---#
    println("--------------------------")
    println("-----Newton 2nd Order-----")
    println("--------------------------")
    println("Converged at iter $N_iter with U = ",U[2])
    PlotSphere(r[:,1,2],r[:,2,2],r[:,3,2],"Newton")
    PlotResiduals(residuals,iterations,"Newton")
    return nothing
end

