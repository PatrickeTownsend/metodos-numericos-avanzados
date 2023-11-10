using LinearAlgebra
using Plots
function PlotSphere(x::Vector,y::Vector,z::Vector, type::String)
    n = 100
    r=0.95
    u = range(-π,π; length=n)
    v = range(0,2π; length=n)
    x_s = r*cos.(u)*sin.(v)'
    y_s = r*sin.(u)*sin.(v)'
    z_s = r*ones(n)*cos.(v)'
    plotly()
    surface(x_s,y_s,z_s, color =:grey, legend=false)
    scatter!(x,y,z, markersize = 3,color=:red)
    if type == "Gradient"
      plt = title!("Gradient descent results")
    elseif  type == "Ipopt"
       plt =title!("Ipopt results")
    elseif type=="Newton"
       plt =title!("Newton 2nd results")
    end
    display(plt)
    #savefig(plt,"Thomsom problem/plots/results_$type.png")
end
function PotentialEnergy(r::Array, N::Int)
    pot = 0.0
    for i = 1:N-1
        for j = i+1:N
            pot += 1/norm(r[i,:] - r[j,:])
        end
    end
    return pot
end

function Gradient(r::Array, i::Int, N::Int)
    grad = zeros(length(r[i,:]))
    for j = 1:N
        if i != j
            grad += (r[j,:] - r[i,:])/norm(r[i,:]-r[j,:])^3
        end
    end
    return grad
end

function hessian(r::Array, i::Int, N::Int)
    H = zeros(length(r[i,:]),length(r[i,:]))
    xᵢ = r[i,1]
    yᵢ = r[i,2]
    zᵢ = r[i,3]
    for j=1:N
        Δx = xᵢ - r[j,1]
        Δy = yᵢ - r[j,2]
        Δz = zᵢ - r[j,3]
        if i!=j
            den = (Δx^2 + Δy^2 + Δz^2)^(5/2)
            H[1,1]+= -(Δy^2 + Δz^2 - 2*Δx^2)/den
            H[2,2]+= -(Δx^2 + Δz^2 - 2*Δy^2)/den
            H[3,3] += -(Δx^2 + Δy^2 - 2*Δz^2)/den
            H[1,2]+= (3*Δx*Δy)/den
            H[1,3]+= (3*Δx*Δz)/den
            H[2,3]+= (3*Δy*Δz)/den
            H[2,1]+= H[1,2]
            H[3,1]+= H[1,3]
            H[3,2]+= H[2,3]
            
        end
    end
    return H
end



function Initialization(N::Int)
    r = 1
    N = N
    x = zeros(N+1)
    y = zeros(N+1)
    z = zeros(N+1)
    Ncount = 1
    a = (4π*r^2)/N
    d = sqrt(a)
    Mϕ =  round(Int, π/d)
    dϕ= π/Mϕ
    dθ = a/dϕ
    for i = 0:Mϕ-1
       ϕ = π*(i+0.5)/Mϕ
       Mθ = round(Int,2π*sin(ϕ)/dθ)
       for j = 0:Mθ-1
           θ = 2π*j/Mθ
           x[Ncount] += r*sin(ϕ)*cos(θ) 
           y[Ncount] += r*sin(ϕ)*sin(θ)
           z[Ncount] += r*cos(ϕ)
           Ncount+=1
        end
    end
    coordinates = hcat(x[1:N],y[1:N],z[1:N])
    return coordinates
end

function InitRandom(N::Int)
    x = 2*rand(N) .- 1
    y = 2*randn(N) .- 1
    z = 2*rand(N) .- 1
    coordinates = hcat(x,y,z)
    return coordinates
end

function PlotResiduals(residuals::Array,N_iterations::Vector,type::String)
    gr()
    plot(N_iterations,[residuals[:,1],residuals[:,2],residuals[:,3] ],
    label=["x residual" "y residual" "z residual"],xlabel="Iterations",ylabel="Residuals",yaxis=:log10)
    plt = Plots.title!("$type residuals")
    savefig(plt,"Thomsom problem/plots/residuals$type.png")
    display(plt)
end