include("optim_ipopt.jl")
include("Functions.jl")
include("gradient_descent.jl")
include("newton_hessian.jl")

#---Set up---#
N = 32 # Number of charges
max_iter = 10000 # Max iteration 
max_tol = 1e-20 # Convergence tolerance
α = 0.01 #Learning rate

#----Initial estimation-----#
r₀_rand = InitRandom(N)
r₀_des = InititDeserno(N)
r₀_fibo = InitFibonacci(N)
r₀_ipopt = optim_ipopt(N,r₀_fibo,max_tol)
InitPotential(r₀_rand,r₀_fibo,r₀_des,r₀_ipopt,N)
init = [r₀_ipopt,r₀_rand,r₀_des,r₀_fibo]
type = ["Ipopt","Random","Deserno","Fibonacci"]
i = 2 #Initial guess selected

#-----Main-------#
@time r,residuals,iterations=gradient_descent(N,init[i],max_tol,max_iter,α,type[i])
#-----Plotting-------#
PlotSphere(r[:,1],r[:,2],r[:,3],type[i])
PlotResiduals(residuals,iterations,type[i])

#-----Write JSON-----#
writeJSON(r,N)
