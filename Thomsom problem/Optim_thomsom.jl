using JuMP
using Ipopt
using Plots
using LinearAlgebra
include("Functions.jl")
model = Model(Ipopt.Optimizer)
set_attribute(model, "print_level",4)
set_optimizer_attribute(model, "max_iter",20000)
set_attribute(model, "print_timing_statistics","yes")

N = 32 #number of charges

@variables(model, begin
    -1 ≤ x[1:N] ≤ 1
    -1 ≤ y[1:N] ≤ 1
    -1 ≤ z[1:N] ≤ 1
end)
coordinates = Initialization(N)
#coordinates = InitRandom(N)
set_start_value.(x,coordinates[:,1])
set_start_value.(y,coordinates[:,2])
set_start_value.(z,coordinates[:,3])
for i = 1:N
    @NLconstraint(model, x[i]^2 + y[i]^2 + z[i]^2 == 1.00000000)
end
@NLobjective(model, Min, sum(sum(((x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2)^(-0.5) for j=i+1:N) for i=1:N-1 ))
optimize!(model)
r = hcat(value.(x),value.(y),value.(z))
U = PotentialEnergy(r,N)
println("----Ipopt----")
println("Potential: $U")
println(r[1,:])
println(r[2,:])
PlotSphere(value.(x),value.(y),value.(z),"Ipopt")
