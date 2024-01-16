module Problems

    export InitialValueProblem
    export IVP1, IVP2


    struct InitialValueProblem
        RHS::Function   # Function returning the right-hand side of the ODE as a column vector
        IC::Vector      # Initial Condition of the ODE. Value is either scalar or column vector.
        tspan::Tuple    # Tuple encoding the integration timespan, i.e. (t0, tf).
    end


    # ----- IVP1: The robertson problem -----

        # Right-hand side of the ODE
        RHS1(t,y::Vector) = [-0.04*y[1] + (1e4)*y[2]*y[3]; 0.04*y[1] - (1e4)*y[2]*y[3] - (3e7)*y[2]^2; (3e7)*y[2]^2];

        # Initial condition of the ODE
        IC1 = [1., 0., 0.];

        # Integration timespan
        tspan1 = (0., 40.);


        # Structure encoding the IVP
        IVP1 = InitialValueProblem(RHS1, IC1, tspan1);


    # ----- IVP2: Prothero-Robinson problem -----

        # Right-hand side of the ODE
        RHS2(t, y) = @. cos(π/4 + t) - (1e6)*(y - sin(π/4 +t))

        # Initial condition of the ODE
        IC2 = [1.];
        
        # Integration timespan
        tspan2 = (0., 1. /10.);
        

        # Structure encoding the IVP
        IVP2 = InitialValueProblem(RHS2, IC2, tspan2);

end