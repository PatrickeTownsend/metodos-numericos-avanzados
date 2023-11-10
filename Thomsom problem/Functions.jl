using LinearAlgebra
using Plots
using JSON
using OrderedCollections
function PlotSphere(x::Vector,y::Vector,z::Vector, type::String)
    n = 100
    r=0.95
    N = length(x)
    u = range(-π,π; length=n)
    v = range(0,2π; length=n)
    x_s = r*cos.(u)*sin.(v)'
    y_s = r*sin.(u)*sin.(v)'
    z_s = r*ones(n)*cos.(v)'
    plotly()
    surface(x_s,y_s,z_s, color =:grey, legend=false)
    scatter!(x,y,z, markersize = 3,color=:red,xlabel="X",ylabel="Y",zlabel="Z")
    plt =title!("Distribution with $type for $N points")
    display(plt)
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


function InititDeserno(N::Int)
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

function writeJSON(r::Array,N::Int)
    data = [OrderedDict{String,Float64}("x" => r[i,1], "y" => r[i,2], "z" => r[i,3]) for i=1:N]
    open("results$N.json","w") do f
       JSON.print(f, data, 4)
    end
end


function InitFibonacci(N::Int)
    r = zeros(N,3)
    Φ = (1+sqrt(5))/2
    for i=0:N-1
        tx = (i+0.5)/(N)
        ty = i/Φ
        θ = acos(2*tx-1)-pi/2
        ϕ = 2*pi*ty
        r[i+1,1] = cos(θ)*cos(ϕ)
        r[i+1,2] = cos(θ)*sin(ϕ)
        r[i+1,3] = sin(θ)
    end
    return r
end

function InitPotential(r_rand::Array,r_fibo::Array,r_des::Array,r_ipopt::Array,N::Int)
    U_rand = PotentialEnergy(r_rand,N)
    U_fibo = PotentialEnergy(r_fibo,N)
    U_des = PotentialEnergy(r_des,N)
    U_ipopt = PotentialEnergy(r_ipopt,N)
    println("---------Initial Potential----------")
    println("Random init   : ",U_rand)
    println("Fibonacci init: ",U_fibo)
    println("Descerno init : ",U_des)
    println("IPOPT init    : ",U_ipopt)
end
