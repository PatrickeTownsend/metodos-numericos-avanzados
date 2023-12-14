using Plots
using LinearAlgebra
f(x) = @. exp(-5*(x-1)^2)
Nx = 10000
xᵢ = LinRange(0,300,Nx)
u_n = zeros(Nx)
u_n1 = zeros(Nx)
Δx = xᵢ[2] - xᵢ[1]
Δt = 0.0001
u_n = f(xᵢ)
for n=1:10000
    for i=2:Nx-1
        u_n1[i] = u_n[i] - (Δt/(2*Δx))*(u_n[i+1] - u_n[i-1])
    end
    u_n = u_n1
end

plot(xᵢ,f(xᵢ))
plot!(xᵢ,u_n1)
xlims!(0,20)
