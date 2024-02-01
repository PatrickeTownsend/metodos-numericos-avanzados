using LinearAlgebra
using OffsetArrays

function W_mult(x_m, x_n, U_m, U_n,k)
    Wt = x_m*U_n[k-1]/k + x_n*U_m[k-1]/k 
    for j=1:(k-1)
       Wt += ((U_m[j-1]/j)*(U_n[k-j-1]/(k-j)))
    end
    return Wt
end

function W_div(x_m, x_n, U_m, U_n, k, Wu)
    Wu = (1/x_n)*(U_m[k-1]/k - sum((U_n[j-1]/j)*Wu[k-j] for j=1:(k)))
    return Wu
end

function W_mult_div(x_m, x_n, U_m, k, Wv)
    Wv = (1/x_m)*(x_n*U_m[k] + sum((U_m[j-1]/j)*U_m[k-j] for j=1:(k)) - sum((U_m[j-1]/j)*Wv[k-j] for j=1:(k)))
    return Wv
end

function initializeX(X,m0,x0, v0)
    X[0,0] = x0[1] #x1
    X[1,0] = x0[2] #x2
    X[2,0] = x0[3] #x3

    X[3,0] = v0[1] #x4
    X[4,0] = v0[2] #x5
    X[5,0] = v0[3] #x6

    X[6,0] = m0 #x7

    X[7,0] = norm(x0)^2 #x8
    X[8,0] = X[7,0]^(3/2) #x9
    X[9,0] = norm(v0)^2 #x10
    X[10,0] = X[9,0]^(1/2) #x11
    return X
end

function initializeU(X,U, μ, T, ṁ)
    U[0,0] = X[3] #u1
    U[1,0] = X[4] #u2
    U[2,0] = X[5] #u3

    U[3,0] = -μ*(X[0]/X[8]) + (T/X[6])*(X[3]/X[10]) #u4
    U[4,0] = -μ*(X[1]/X[8]) + (T/X[6])*(X[4]/X[10]) #u5
    U[5,0] = -μ*(X[2]/X[8]) + (T/X[6])*(X[5]/X[10]) #u6

    U[6,0] = ṁ #u7

    U[7,0] = 2*X[0]*X[3] + 2*X[1]*X[4] + 2*X[2]*X[5] #u8
    U[8,0] = (3/2)*(X[8]*U[7,0])/(X[7]) #u9

    U[9,0] = 2*X[3]*U[3,0] + 2*X[4]*U[4,0] + 2*X[5]*U[5,0] #u10
    U[10,0] = (1/2)*(X[10]*U[9,0])/(X[9]) #u11

    return U
end

function initializeW(X,U,W)

    W[0,0] = X[0]/X[8] # W4,1
    W[1,0] = X[3]/X[10] # W4,2

    W[2,0] = X[1]/X[8] # W5,1
    W[3,0] = X[4]/X[10] # W5,2

    W[4,0] = X[2]/X[8] # W6,1
    W[5,0] = X[5]/X[10] # W6,2

    W[6,0] = X[0]*U[0] # W8,1
    W[7,0] = X[1]*U[1] # W8,2
    W[8,0] = X[2]*U[2] # W8,3

    W[9,0] = X[8]*U[7]/X[7] # W9

    W[10,0] = X[3]*U[3] # W10,1
    W[11,0] = X[4]*U[4] # W10,2
    W[12,0] = X[5]*U[5] # W10,3

    W[13,0] = X[10]*U[9]/X[9] # W11

    return W
end

function Recursive(X,U,W,k_max,j)
    for k = 1:k_max-1
        # Calculate W
        W[0,k] = W_div(X[0,j],X[8,j],U[0,:],U[8,:],k, W[0,:]) # W4,1
        W[1,k] = W_div(X[3,j],X[10,j],U[3,:],U[10,:],k, W[1,:]) # W4,2

        W[2,k] = W_div(X[1,j],X[8,j],U[1,:],U[8,:],k, W[2,:]) # W5,1
        W[3,k] = W_div(X[4,j],X[10,j],U[3,:],U[10,:],k, W[3,:]) # W5,2

        W[4,k] = W_div(X[2,j],X[8,j],U[3,:],U[8,:],k, W[4,:]) # W6,1
        W[5,k] = W_div(X[5,j],X[10,j],U[5,:],U[10,:],k, W[5,:]) # W6,2

        W[6,k] = W_mult(X[0,j],X[3,j],U[0,:],U[3,:],k) # W8,1
        W[7,k] = W_mult(X[1,j],X[4,j],U[1,:],U[4,:],k) # W8,2
        W[8,k] = W_mult(X[2,j],X[5,j],U[2,:],U[5,:],k) # W8,3

        W[9,k] = W_mult_div(X[7,j],X[8,j],U[7,:],k, W[9,:]) # W9

        W[10,k] = (X[3,j]*U[3,k] + (U[3,k-1]/k)*U[3,k])
        W[11,k] = (X[4,j]*U[4,k] + (U[4,k-1]/k)*U[4,k])
        W[12,k] = (X[5,j]*U[5,k] + (U[5,k-1]/k)*U[5,k])

        for j = 1:(k-1)
           W[10,k] += ((U[3,j-1]/j)*U[3,k-j]) # W10,1
           W[11,k] += ((U[4,j-1]/j)*U[4,k-j]) # W10,2
           W[12,k] += ((U[5,j-1]/j)*U[5,k-j]) # W10,3
        end
        W[13,k] = W_mult_div(X[9,j],X[10,j],U[9,:],k, W[13,:]) # W11

        # Calculate U(k)
        U[10,k] = 0.5*W[13,k] # U11
        U[9,k] = 2*W[10,k] + 2*W[11,k] + 2*W[10,k] # U10
        U[8,k] = (3/2)*W[9,k] # U9
        U[7,k] = 2*W[6,k] + 2*W[7,k] + 2*W[8,k] # U8
        U[6,k] = 0 # U7
        U[5,k] = -μ*W[4,k] + (T/X[6,j])*W[5,k] # U6
        U[4,k] = -μ*W[2,k] + (T/X[6,j])*W[3,k] # U5
        U[3,k] = -μ*W[0,k] + (T/X[6,j])*W[1,k] # U4
        U[2,k] = U[5,k-1]/k
        U[1,k] = U[4,k-1]/k
        U[0,k] = U[3,k-1]/k
    end
    return U, W
end













