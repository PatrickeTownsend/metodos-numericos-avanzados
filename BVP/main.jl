using LinearAlgebra
using Plots

# Set problem parameters:
k  = 1.     # k: square of oscillator's frequency
y0 = 1.     # y0: initial value for y(t0)


# Timespan
t0 = 0.
tf = 10.


# Time discretization
N = 100
h = @. (tf-t0) / N


# Initialize matrix A and vector b
A = @. ones(Float64, N, N) * 0.
b = @. ones(Float64, N, 1) * 0.


# Fill in vector b
b[1] = - y0


# Fill in matrix A:

    # Diagonal bands
    for i in 2:(N-1)
        A[i,i-1] =  1.
        A[i,i]   = @. -2 + k*h^2
        A[i,i+1] =  1.
    end

    # i = 1 (Impose boundary value at t1 = t0 + h)
    A[1,1] = @. k*h^2 - 2
    A[1,2] = 1.

    # i = N (Impose boundary value at tf)
    A[N,N-1] = -h
    A[N,N]   =  h


# Solve Linear System: Get solution y(t) at discrete instants of time
y = A\b


# Add value at t0 to the solution
y = [y0; y]


# Get vector for equi-spaced, discretized time values
t = ones(Float64, N, 1)
for i in 1:N
    t[i] = i .* h
end
t = [t0; t]


# Plot solution
p = plot(t, y)
display(p)
plot!(t, y, seriestype=:scatter, ms=2)

println("Done!")