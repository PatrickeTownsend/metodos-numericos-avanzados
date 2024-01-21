module Problems

    export InitialValueProblem
    export IVP1, IVP2


    struct InitialValueProblem
        RHS::Function   # Function returning the right-hand side of the ODE as a column vector
        IC::Vector      # Initial Condition of the ODE. Value is either scalar or column vector.
        tspan::Tuple    # Tuple encoding the integration timespan, i.e. (t0, tf).
        Jacobian::Function
    end


    # ----- IVP1: The robertson problem -----

        # Right-hand side of the ODE
        RHS1(t,y::Vector) = [-0.04*y[1] + (1e4)*y[2]*y[3]; 0.04*y[1] - (1e4)*y[2]*y[3] - (3e7)*y[2]^2; (3e7)*y[2]^2];

        # Initial condition of the ODE
        IC1 = [1., 0., 0.];

        # Integration timespan
        tspan1 = (0., 40.);

        function J1(y)
            dF1dy1 = (-0.04)
            dF1dy2 =  (1e4)*y[3]
            dF1dy3 = (1e4)*y[2]
    
            dF2dy1 = (0.04)
            dF2dy2 = ((-1e4)*y[3] - 2*(3e7)*y[2])
            dF2dy3 = (-1e4)*y[2]
    
            dF3dy1 = 0
            dF3dy2 = (2*3e7)*y[2]
            dF3dy3 = 0
    
            J = [dF1dy1 dF1dy2 dF1dy3;
                 dF2dy1 dF2dy2 dF2dy3;
                 dF3dy1 dF3dy2 dF3dy3]
    
            return J
        end

       


        # Structure encoding the IVP
        IVP1 = InitialValueProblem(RHS1, IC1, tspan1, J1);


    # ----- IVP2: Prothero-Robinson problem -----

        # Right-hand side of the ODE
        RHS2(t, y) = @. cos(π/4 + t) - (1e6)*(y - sin(π/4 +t))

        # Initial condition of the ODE
        IC2 = [1.];
        
        # Integration timespan
        tspan2 = (0., 1. /10.);

        J2(y) = [-1e6]

        # Structure encoding the IVP
        IVP2 = InitialValueProblem(RHS2, IC2, tspan2, J2);

end