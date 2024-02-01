# ----- Initialization -----

    # Load need libraries
    using Plots
    using LinearAlgebra
    
    # Include externally declared functions
    include("Problems.jl")
    include("RungeKutta.jl")

    using .Problems
    using .RungeKutta



# --- Integrate given problem with an Explicit Runge-Kutta method ---

    #Set the IVP to be solved
    IVP = Problems.IVP2;
    
    # Set up the integration parameters:
    order =  3;    # Set the order of the Explicit Runge-Kutta integrator
    h =  0.01;      # Set the stepsize

    # Perform integration
    t, y, f = RungeKutta.Explicit(IVP, order, h)
    display(y[:,end])

    p1 = plot(t', y', title="Numerical Solution (h=$h)", label="RK$order")
    xlabel!("t")
    ylabel!("y")
    display(p1)

    # The integrator returns the computed numerical solution stored in the following variables:
    #  --  t = [t₀ t₁ t₂ t₃ ...]    where tᵢ are scalars representing time instants where the solution is provided
    #  --  y = [y₀ y₁ y₂ y₃ ...]    where yᵢ are column vectors representing the solution y(tᵢ) at times tᵢ
    #  --  f = [f₀ f₁ f₂ f₃ ...]    where fᵢ are column vectors representing the derivatives f(tᵢ, yᵢ)


#     # Compute Analytical Solution
#     yexact = IVP.Sol.(t)


#     # We can also integrate with a different order methods, for the sake of comparison
#     t1, y1, f1 = Adams.Bashforth(IVP, 1, h)
#     t2, y2, f2 = Adams.Bashforth(IVP, 2, h)
#     t3, y3, f3 = Adams.Bashforth(IVP, 3, h)
#     t4, y4, f4 = Adams.Bashforth(IVP, 4, h)


# # --- Plot results ---

#     println("Plotting...")

#     # Plot provided-order solution
#     p1 = plot(t', y', title="Numerical Solution (h=$h)", label="RK$order")
#     xlabel!("t")
#     ylabel!("y")

#     # Compute errors wrt analytic solution
#     function RelativeError(IVP, t, y, yexact)

#         if length(IVP.IC) == 1
#             e = abs.((y .- yexact) ./ yexact);
#         else
#             e = 0. .* t;
#             for i = 1:length(t)
#                 e[i] = norm( (y[:,i] .- yexact[i]) ./ yexact[i] );
#             end
#         end

#         return e
#     end

#     e1 = RelativeError(IVP, t1,  y1,  yexact);
#     e2 = RelativeError(IVP, t2, y2, yexact);
#     e3 = RelativeError(IVP, t3, y3, yexact);
#     e4 = RelativeError(IVP, t4, y4, yexact);
    
#     # Plot numerical error vs the analytical solution
#     p2 = plot(t1', e1', yaxis=:log10, title="Relative Error of the Numerical Solution (h=$h)", label="RK1")
#     plot!(t2', e2', label="RK2")
#     plot!(t3', e3', label="RK3")
#     plot!(t4', e4', label="RK4")
#     xlabel!("t")
#     ylabel!("error")

#     p = plot(p1, p2, layout=(2, 1))
#     display(p)

    println("Done!\n")

    
