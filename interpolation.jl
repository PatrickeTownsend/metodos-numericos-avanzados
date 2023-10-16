using Plots
N = 100 #Numero de puntos
x_i = LinRange(-10,10,N)
f(x) = 1 ./(1 .+25 .*x.^2)
fx = f(x_i)
function halfPoint(x_half, x, fx)
    for i in eachindex(x)
        if i == N
            lower = x[i]
            upper = x[N]
            y_u = fx[N]
            y_l = fx[i]
        else
            lower = x[i]
            upper = x[i+1]
            y_u = fx[i+1]
            y_l = fx[i]
        end
        if lower ≤ x_half ≤ upper
            m = (y_u - y_l)/(upper - lower)
            return y = m*(x_half - upper) + y_u
        end
    end
end
x_need = 0.3
y = halfPoint(x_need, x_i, fx)
exact = f(x_need)
error = abs(exact-y)
println("valor objetivo $x_need")
println("error = $error")
println("valor interpolado = $y")
println("valor exacto = $exact")
#scatter(x_i,f(x_i))