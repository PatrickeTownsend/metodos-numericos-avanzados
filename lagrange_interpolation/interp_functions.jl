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


function ChebyshevNode(x_init, x_end, n)
    x_ch = zeros(n)
    for i =1:n
        x_ch[i] = ((x_end-x_init)/2)*cos((2*i-1)*π/(2*n)) + (x_init + x_end)/2
    end
    return x_ch
end

function Plot_func(xₗ, yₗ, Pᵢ)
    Plots.plot(xₗ , yₗ, label="f(x) exacta")
    plot= Plots.scatter!(xₗ , Pᵢ, label="F(x) numerica")
    Plots.xlabel!("xₗ")
    ylabel!("y")
    title!("Interpolacion")
    display(plot)
end

function distribucionB(x_init, x_end, n)
    x_ch = zeros(n+1)
    x_ch[1] += cos(π*((2*0+1)/(2*n+2))) 
    for i =1:n
        x_ch[i+1] = cos(π*((2*i+1)/(2*n+2))) 
    end
    return x_ch
end

function  distribucionC(x_init,x_end,n)
    x_ch = zeros(n+1)
    x_ch[1]+=cos(π*0/n)
    for i =1:n
        x_ch[i+1] = cos(π*i/n) 
    end
    return x_ch
end
