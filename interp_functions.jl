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
        poly = poly + Li*yᵢ
    end
    return poly,Li
end

function epsilon(x,xᵢ)
    N = length(xᵢ)
    π = prod(x - tᵢ for tᵢ in xᵢ)/(factorial(big(N+1)))
    return π
end