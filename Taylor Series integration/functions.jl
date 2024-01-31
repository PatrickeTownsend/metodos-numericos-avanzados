using LinearAlgebra
using OffsetArrays

function W_mult(x_m, x_n, U_m, U_n,k)
    Wt = x_m*U_n[k-1]/k + x_n*U_m[k-1]/k + sum((U_m[j-1]/j)*(U_n[k-j-1]/(k-j)) for j=1:(k-1))
    return Wt
end

function W_div(x_m, x_n, U_m, U_n,k)
    Wu = (1/x_n)*(U_m[k-1]/k - sum((U_n[j-1]/j)*Wu[k-j] for j=1:(k)))
    return Wu
end

function W_mult_div(x_m, x_n, U_m, U_n,k)
    Wv = (1/x_m)*(x_n*U_m[k] + sum((U_m[j-1]/j)*U_m[k-1] for j=1:(k)) - sum((U_m[j-1]/j)*Wv[k-j] for j=1:(k)))
    return Wv
end

function initializeX(X,m0,x0, v0)
    X[0,0] += x0[1] #x1
    X[1,0] += x0[2] #x2
    X[2,0] += x0[3] #x3

    X[3,0] += v0[1] #x4
    X[4,0] += v0[2] #x5
    X[5,0] += v0[3] #x6

    X[6,0] += m0 #x7

    X[7,0] += norm(x0)^2 #x8
    X[8,0] += X[7,0]^(3/2) #x9
    X[9,0] += norm(v0)^2 #x10
    X[10,0] += X[9,0]^(1/2) #x11
    return X
end

function initializeU(X,U, μ, T, ṁ)
    U[0,0] += X[3,0] #u1
    U[1,0] += X[4,0] #u2
    U[2,0] += X[5,0] #u3

    U[3,0] += -μ*(X[0,0]/X[8,0]) + (T/X[6,0])*(X[3,0]/X[10,0]) #u4
    U[4,0] += -μ*(X[1,0]/X[8,0]) + (T/X[6,0])*(X[4,0]/X[10,0]) #u5
    U[5,0] += -μ*(X[2,0]/X[8,0]) + (T/X[6,0])*(X[5,0]/X[10,0]) #u6

    U[6,0] += ṁ #u7

    U[7,0] += 2*X[0,0]*X[3,0] + 2*X[1,0]*X[4,0] + 2*X[2,0]*X[5,0] #u8
    U[8,0] += (3/2)*(X[8,0]*U[7,0])/(X[7,0]) #u9

    U[9,0] += 2*X[3,0]*U[3,0] + 2*X[4,0]*U[4,0] + 2*X[5,0]*U[5,0] #u10
    U[10,0] += (1/2)*(X[10,0]*U[9,0])/(X[9,0]) #u11

    return U
end

function initializeW(X,U,W)

    W[0,0] += X[0,0]/X[8,0] # W4,1
    W[1,0] += X[3,0]/X[10,0] # W4,2

    W[2,0] += X[1,0]/X[8,0] # W5,1
    W[3,0] += X[4,0]/X[10,0] # W5,2

    W[4,0] += X[2,0]/X[8,0] # W6,1
    W[5,0] += X[5,0]/X[10,0] # W6,2

    W[6,0] += X[0,0]*U[0,0] # W8,1
    W[7,0] += X[1,0]*U[1,0] # W8,2
    W[8,0] += X[2,0]*U[2,0] # W8,3

    W[9,0] += X[8,0]*U[7,0]/X[7,0] # W9

    W[10,0] += X[3,0]*U[3,0] # W10,1
    W[11,0] += X[4,0]*U[4,0] # W10,2
    W[12,0] += X[5,0]*U[5,0] # W10,3

    W[13,0] += X[10,0]*U[10,0]/X[9,0] # W11

    return W
end

function Recursive(X,U,W,k_max,i)
    for k = 1:k_max
        W[0,k] += W_div(X[0,j],X[8,j],U[0,:],U[8,:],k) # W4,1
        W[1,k] += W_div(X[3,j],X[10,j],U[3,:],U[10,:],k) # W4,2

        W[2,k] += W_div(X[1,j],X[8,j],U[1,:],U[8,:],k) # W5,1
        W[3,k] += W_div(X[4,j],X[10,j],U[3,:],U[10,:],k) # W5,2

        W[4,k] += W_div(X[2,j],X[8,j],U[3,:],U[8,:],k) # W6,1
        W[5,k] += W_div(X[5,j],X[10,j],U[5,:],U[10,:],k) # W5,1

        





