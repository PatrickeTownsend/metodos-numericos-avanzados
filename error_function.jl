using Plots
N = 50 #Numero de puntos
xi = -10
xf = 10
Δx = (xf-xi)/(N*100)
x_i = LinRange(-10,10,N)
f(x) = 1 ./(1 .+25 .*x.^2)
π = zeros(N)
for i = 1:N
    x = x_i[i]
    for j = 0:N
       π[i] *= ((x/Δx) - j)
    end
end

ϵ = @. ((Δx^(N+1))/factorial(big(N+1)))*π
plot(x_i, π)
println(π)