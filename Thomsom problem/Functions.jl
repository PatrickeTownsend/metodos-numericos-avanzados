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
