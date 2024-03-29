using LinearAlgebra
using Plots
using JSON
using OrderedCollections
"""
Plots the sphere with the charges equidistributed

# INPUT
- `x::Vector` X coordinates of the charges
- `y::Vector` Y coordinates of the charges
- `z::Vector` Z coordinates of the charges
- `type::String` type of distribution

"""
function PlotSphere(x::Vector,y::Vector,z::Vector, type::String)
    n = 100
    r=0.95
    N = length(x)
    u = range(-π,π; length=n)
    v = range(0,2π; length=n)
    x_s = r*cos.(u)*sin.(v)'
    y_s = r*sin.(u)*sin.(v)'
    z_s = r*ones(n)*cos.(v)'
    plotlyjs()
    surface(x_s,y_s,z_s, color =:grey, legend=false)
    scatter!(x,y,z, markersize = 4,color=:red,xlabel="X",ylabel="Y",zlabel="Z")
    plt =title!("$type distribution for $N points")
    #savefig(plt,"Thomsom problem/plots/$type Distribution $N.png")
    display(plt)
end
"""
Calculates the potential energy of the whole configuration

# INPUT
- `r::Array` Array with coordinates of each charges
- `N::Int` Number of charges to place

# OUTPUT 
- `pot::Float64` Potential energy
"""
function PotentialEnergy(r::Array, N::Int)
    pot = 0.0
    for i = 1:N-1
        for j = i+1:N
            pot += 1/norm(r[i,:] - r[j,:])
        end
    end
    return pot
end
"""
Calculates the gradient of the function for a charge i
# INPUT 
- `r::Array` Array with coordinates of each charges
- `i::Int` Number of the charge
- `N::Int` Number of charges to place

# OUTPUT 
- `grad::Vector` Gradiente vector of charge i 
"""
function Gradient(r::Array, i::Int, N::Int)
    grad = zeros(length(r[i,:]))
    for j = 1:N
        if i != j
            grad += (r[j,:] - r[i,:])/norm(r[i,:]-r[j,:])^3
        end
    end
    return grad
end
"""
Deserno initial configuration

# INPUT
- `N::Int` Number of charges

# OUTPUT 
- `coordinates::Array` Array with the initial configuration 
"""

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
"""
Generate a random initial configuration
# INPUT
- `N::Int` Number of charges

# OUTPUT 
- `coordinates::Array` Array with the initial configuration 
"""
function InitRandom(N::Int)
    θ = rand(N)*2pi 
    ϕ = acos.(2 .*rand(N) .- 1)
    x = @. sin(ϕ)*cos(θ)
    y = @. sin(ϕ)*sin(θ)
    z = @. cos(ϕ)
    coordinates = hcat(x,y,z)
    return coordinates
end
"""
Plots the residuals and derivatives of the iterations
"""
function PlotResiduals(derivatives::Array,N_iterations::Vector,grad::Vector, type::String,init::String)
    gr()
    plot(N_iterations,[derivatives[:,1],derivatives[:,2],derivatives[:,3] ],
    label=["x residual" "y residual" "z residual"],xlabel="Iterations",ylabel="Derivatives",yaxis=:log10,
    linewidth=2)
    plt = Plots.title!("$type evolution $init initialization")
    plt2 = plot(N_iterations,grad, xlabel="Iterations",ylabel="Residual",title="Iteration residuals $init initialization",
    legend=false)
    #savefig(plt,"Thomsom problem/plots/derivatives$type $init.png")
    #savefig(plt2,"Thomsom problem/plots/residual$type $init.png")
    display(plt)
    display(plt2)
end
"""
Writes a JSON file with each charge position vector
"""
function writeJSON(r::Array,N::Int)
    data = [OrderedDict{String,Float64}("x" => r[i,1], "y" => r[i,2], "z" => r[i,3]) for i=1:N]
    open("Thomsom problem/results$N.json","w") do f
       JSON.print(f, data, 4)
    end
end
"""
Generate a Fibonacci Lattice configuration
# INPUT
- `N::Int` Number of charges

# OUTPUT 
- `r::Array` Array with the initial configuration 
"""

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
"""
Print the potential for each initial distribution

# INPUT 
- `r_rand::Array` Random configuration coordinates
- `r_fibo::Array` Fibonacci configuration coordinates
- `r_des::Array` Deserno configuration coordinates
- `r_ipopt::Array` IPOPT configuration coordinates
- `N::Int` Number of charges
"""
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
