using LinearAlgebra

function F(x) 
    F = zeros((length(x)))
    F[1] = -x[1]^3 + x[2]
    F[2] = x[1]^2 + x[2]^2 -1
    return F
end
function J(x)
    [
        -3*x[1]^2 1;
        2*x[1] 2*x[2]
    ]
end
x = [1,2]
Func = F(x)
Jac = J(x)
Tolx = 1.2e-12
TolF = 1.2e-12
for i = 1:20
    δx = -Jac \ Func
    x += δx
    
    println("Iter $i: $x")
    if (norm(δx) ≤ Tolx) && (norm(Func)≤TolF)
        break
    end
    Func = F(x)
    Jac = J(x)
end

