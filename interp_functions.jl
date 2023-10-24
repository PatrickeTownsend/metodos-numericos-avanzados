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

function ChebyshevNode(n)
    x_ch = zeros(n)
    for i =1:n
        x_ch[i] = cos((2*i-1)*π/(2*n))
    end
    return x_ch
end

function Plot_func(xₗ, yₗ, Pᵢ)
    plot(xₗ , yₗ, label="f(x) exacta")
    scatter!(xₗ , Pᵢ, label="F(x) numerica")
    xlabel!("xₗ")
    ylabel!("y")
    title!("Interpolacion")
end