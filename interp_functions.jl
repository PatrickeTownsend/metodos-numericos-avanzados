function Lagrange(x ,xᵢ ,yᵢ)
    N = length(xᵢ)
    poly = 0.0
    for i = 1:N
        Li = 1.0
        for j = 1:N
            if j != i
                Li *= (x - xᵢ[j])/(xᵢ[i] - xᵢ[j])
            end
        end
        poly = poly + Li*yᵢ[i]
    end
    return poly
end

function epsilon(x,xᵢ)
    N = length(xᵢ)
    π = prod(x - tᵢ for tᵢ in xᵢ)/(factorial(big(N+1)))
    return π
end

function epsilon_eq(x,xᵢ)
    N = length(xᵢ)
    Δx = (xᵢ[N]-xᵢ[1])/N
    cte = (Δx^(N+1))/(factorial(big(N+1)))
    π = prod(x/Δx - xᵢ[j] for j =1:N)
    return cte*π
end

function ChebyshevNode(x_init, x_end, n)
    x_ch = zeros(n)
    for i =1:n
        x_ch[i] = ((x_end-x_init)/2)*cos((2*i-1)*π/(2*n)) + (x_init + x_end)/2
    end
    return x_ch
end

function Plot_func(xₗ, yₗ, Pᵢ)
    Plots.plot(xₗ , yₗ, label="f(x) exacta")
    Plots.scatter!(xₗ , Pᵢ, label="F(x) numerica")
    Plots.xlabel!("xₗ")
    Plots.ylabel!("y")
    Plots.title!("Interpolacion")
end

#function LagrangeBary(x ,xᵢ ,yᵢ)
