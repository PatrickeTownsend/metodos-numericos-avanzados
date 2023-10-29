function PlotSphere(x,y,z)
    n = 100
    u = range(-π,π; length=n)
    v = range(0,2π; length=n)
    x_s = cos.(u)*sin.(v)'
    y_s = sin.(u)*sin.(v)'
    z_s = ones(n)*cos.(v)'
    plotly()
    surface(x_s,y_s,z_s)
    scatter!(x,y,z)
    title!("Gradient descent results")
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
            grad += (r[j,:] - r[i,:])/norm(r[i,:]-r[j,:])
        end
    end
    return grad
end

function Initialization(N::Int)
    r = 1
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
    coordinates = hcat(x[1:Ncount-1],y[1:Ncount-1],z[1:Ncount-1])
    return coordinates
end