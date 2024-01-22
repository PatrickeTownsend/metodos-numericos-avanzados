using LinearAlgebra
using Plots
function Regression(f::Function, type::String, p::Int, x::LinRange, ω::Float64)
    y = f(x)
    N = size(x)[1]

    if type == "polynomial"
        X = @. ones(Float64, N, 1)
        for n = 1:p
            X = @. [X x^n]
        end
        A = (X'*X) \ (X'*y)
        f_reg = 0

        for i=1:p+1
            f_reg = @. f_reg + A[i].*x.^(i-1)
        end

    elseif type == "fourier"
        X = @. ones(Float64, N, 1)
        for i = 1:p
            X = @. [X cos(i*x*ω) sin(i*x*ω)]
        end
        A = (X'*X) \ (X'*y)
        f_reg = A[1]

        for i=1:p
            f_reg = @. f_reg + A[2*i]*cos(i*ω*x) + A[2*i+1]*sin(i*ω*x)
        end
    end
        

    return f_reg
end

a = 0
b = 100
N = 1000
f(x) = @. 1 + 0.5*x + 5. *randn() + 10*sin(0.3*x) - exp(0.045*x)
x = LinRange(a,b,N)
f_reg = Regression(f,"polynomial", 4, x, 0.3)
f_fourier = Regression(f,"fourier", 150, x, 0.01)

plot(x, f(x), seriestype=:scatter, ms=2, label="Data")
plot!(x,f_reg, linewidth=4, label="regression")
plot!(x, f_fourier, lw=4, label="fourier")