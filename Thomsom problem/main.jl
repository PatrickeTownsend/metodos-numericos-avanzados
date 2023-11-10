include("optim_ipopt.jl")
include("Functions.jl")
include("gradient_descent.jl")
include("newton_hessian.jl")
#---Set up---#
N = 32 # Number of charges
max_iter = 10000 # Max iteration 
max_tol = 1e-20 # Convergence tolerance
α = 0.01 #Learning rate
r₀ = InitRandom(N)
#r₀ = Initialization(N)
@time optim_ipopt(N,r₀,max_tol)
@time gradient_descent(N,r₀,max_tol,max_iter,α)
#newton2nd(N,r₀,max_tol,max_iter,α)