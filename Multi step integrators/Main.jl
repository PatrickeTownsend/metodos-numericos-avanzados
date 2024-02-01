# ----- Initialization -----

    # Load need libraries
    using Plots
    #using ProgressBars
    
    # Include externally declared functions
    include("Problems.jl")
    include("Adams.jl")

    import .Problems
    import .Adams



# --- Integrate given problem with an Adams method ---

    #Set the IVP to be solved
    IVP = Problems.IVP1;
    
    # Set up the integration parameters:
    k =  4;         # Set the number of steps that the Adams integrator employs at each integration step
    h =  0.01;      # Set the stepsize

    # Perform integration
    t, y, f = Adams.Bashforth(IVP, k, h)

    # The integrator returns the computed numerical solution stored in the following variables:
    #  --  t = [t₀ t₁ t₂ t₃ ...]    where tᵢ are scalars representing time instants where the solution is provided
    #  --  y = [y₀ y₁ y₂ y₃ ...]    where yᵢ are column vectors representing the solution y(tᵢ) at times tᵢ
    #  --  f = [f₀ f₁ f₂ f₃ ...]    where fᵢ are column vectors representing the derivatives f(tᵢ, yᵢ)


    # Compute Analytical Solution
    println("Computing analytic  solution...")
    yexact = IVP.Sol.(t)


    # We can also integrate with a higher-order methods, for the sake of comparison
    k2 = k + 1
    t2, y2, f2 = Adams.Bashforth(IVP, k2, h)

    k3 = k + 2
    t3, y3, f3 = Adams.Bashforth(IVP, k3, h)

    k4 = k + 3
    t4, y4, f4 = Adams.Bashforth(IVP, k4, h)


# --- Plot results ---

    println("Plotting...")

    # Plot lowest-order solution
    p1 = plot(t', y', title="Numerical Solution", label="AB$k")
    xlabel!("t")
    ylabel!("y")

    # Plot highest-order solution
    p2 = plot(t', y4', title="Numerical Solution", label="AB$k4")
    xlabel!("t")
    ylabel!("y")

    # Compute errors wrt analytic solution
    function RelativeError(IVP, t, y, yexact)

        if length(IVP.IC) == 1
            e = abs.((y .- yexact) ./ yexact);
        else
            e = 0. .* t;
            for i = 1:length(t)
                e[i] = sqrt( sum( abs2.(y[:,i]) )) ;
            end
        end

        return e
    end

    e  = RelativeError(IVP, t,  y,  yexact);
    e2 = RelativeError(IVP, t2, y2, yexact);
    e3 = RelativeError(IVP, t3, y3, yexact);
    e4 = RelativeError(IVP, t4, y4, yexact);
    
    # Plot numerical error vs the analytical solution
    p3 = plot(t', e', title="Relative Error of the Numerical Solution", label="AB$k")
    plot!(t2', e2', label="AB$k2")
    plot!(t3', e3', label="AB$k3")
    plot!(t4', e4', label="AB$k4")
    xlabel!("t")
    ylabel!("error")

    p = plot(p1, p2, p3, layout=(3, 1))
    display(p)

    println("Done!\n")

