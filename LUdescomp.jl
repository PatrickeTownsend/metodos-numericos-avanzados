using LinearAlgebra

function LUfactor(A)
    N = size(A)[1]
    L = Matrix{Float64}(I,N,N)
    U = zeros(N,N)

    for j = 1:N
       for i = 1:j
           U[i,j] = A[i,j] - sum(L[i,k]*U[k,j] for k=1:i)
       end
       for i = (j):(N)
           L[i,j] = (1/U[j,j])*(A[i,j]- sum(L[i,k]*U[k,j] for k=1:j))
       end
       L[j,j] = 1.0
    end
    return L, U
end

function LU_linear_solve(L,U,b)
    N = size(L)[1]
    y = zeros(N)
    x = zeros(N)
    y[1] = b[1]/L[1,1]
    for i=2:N
        y[i] = (1/L[i,i])*(b[i] - sum(L[i,j]*y[j] for j=1:i))
    end
    x[N] = y[N]/U[N,N]
    for i = (N-1):-1:1
        x[i] = (1/U[i,i])*(y[i] - sum(U[i,j]*x[j] for j=i:N))
    end
    return x
end


