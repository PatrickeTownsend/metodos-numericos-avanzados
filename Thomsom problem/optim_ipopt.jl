using JuMP
using Ipopt
include("Functions.jl")

function optim_ipopt(N,r₀,max_tol)
    model = Model(Ipopt.Optimizer)
    set_attribute(model, "print_level",0)
    set_optimizer_attribute(model, "max_iter",20000)
    set_optimizer_attribute(model, "acceptable_tol",max_tol)
    @variables(model, begin
       -1 ≤ x[1:N] ≤ 1
       -1 ≤ y[1:N] ≤ 1
       -1 ≤ z[1:N] ≤ 1
    end)
    set_start_value.(x,r₀[:,1])
    set_start_value.(y,r₀[:,2])
    set_start_value.(z,r₀[:,3])
    for i = 1:N
        @NLconstraint(model, x[i]^2 + y[i]^2 + z[i]^2 == 1.00000000)
    end
    @NLobjective(model, Min, sum(sum(((x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2)^(-0.5) for j=i+1:N) for i=1:N-1 ))
    optimize!(model)
    r = hcat(value.(x),value.(y),value.(z))
    return r
end
