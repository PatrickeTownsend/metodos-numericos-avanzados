using Plots
x = range(0,100,length=1000)
y = @. 1 + 0.5 * x + 5. * randn() + 10 *sin(0.3 *x) - exp(0.045*x) #@. es la amcro para poner el punto en todos los operadores

plot(x,y, seriestype=:scatter, ms=2, label="data")


X = [ones(Float64, 1000,1) x x.^2 x.^3 x.^4]
A = (X'*X) \ (X'*y)

yr =  @. A[1] + A[2]*x + A[3]*x^2 + A[4]*x^3 + A[5]*x^4

plot!(x,yr, linewidth=4, label="regresion")

#-----Regresion truncated fourier-----#

# a0 + a1*cos(x) + b1*sin(x) + a2*cos(2x) + b2*sin(2x) + ...
ω = 0.01
X = @. ones(Float64, 1000,1)
for i = 1:150
    X = @. [X cos(i*x*ω) sin(i*x*ω)]
end
A = (X'*X) \ (X'*y)

yr = A[1]
for i = 1:150
    yr = @. yr + A[2*i]*cos(i*ω*x) + A[2*i+1]*sin(i*ω*x)
end

plot!(x,yr, lw=4, label="fourier")