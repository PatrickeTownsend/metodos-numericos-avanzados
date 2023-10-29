using Plots
using LinearAlgebra
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
end
N = 32
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
println(Ncount)
PlotSphere(x[1:Ncount-1],y[1:Ncount-1],z[1:Ncount-1])
coordinates = hcat(x,y,z)


