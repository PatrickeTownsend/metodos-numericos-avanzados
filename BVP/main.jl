using LinearAlgebra
using Plots
t0 = 0.
tf = 10.
N = 100
k = 1.
h = (tf-t0)/(N+1)

A = zeros(Float64, N, N)
B = zeros(Float64, N, 1)
B[1] = -h - k*h^2

#bands
for i=2:(N-1)
    A[i,i-1] = 1
    A[i,i] = -2. + k*h^2
    A[i,i+1] = 1.
end

# i=1
A[1,1] =  k*h^2 - 2*h
A[1,2] = 1.

# i = N
A[N,N-1] = -h
A[N,N] = h

y = A\B
t = zeros(Float64,N,1)
for i in eachindex(t)
    t[i] = i .* h
end

plot(t,y)


