using Plots
#--- Definicion del intervalo---#
xi = -10
xf = 10
#---Numero de puntos---#
N = 20 # Puntos a interpolar
n = 1000 # Puntos para graficar
x_l = LinRange(xi,xf,n)
x_i = LinRange(xi,xf,N)
fi  = zeros(n)
poly = 0
f(x) = @.  1/(1 +25 *x^2)
for k in eachindex(x_l)
    x = x_l[k]
    poly = 0.0
    for i = 1:N
        Li = 1.0
        for j = 1:N
            if j != i
                Li *= (x - x_i[j])/(x_i[i] - x_i[j])
            end
        end
        poly = poly + Li*f(x_i[i])
    end
    fi[k] += poly
end
f_exact = f(x_l)
error = abs.(f_exact .- fi)
plot(x_l,f_exact)
scatter!(x_l,fi)
scatter(x_l, log.(error))

