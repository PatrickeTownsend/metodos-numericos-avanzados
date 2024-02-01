# This file defines multiple "Initial Value Problems" (IVP).
# Thus, the file acts as a library that encodes a series of
# prescribed IVPs that can later be used with different
# time steps and integrators.

module Problems

    export InitialValueProblem
    export IVP1, IVP2


    struct InitialValueProblem
        RHS::Function   # Function returning the right-hand side of the ODE as a column vector
        IC::Vector      # Initial Condition of the ODE. Value is either scalar or column vector.
        tspan::Tuple    # Tuple encoding the integration timespan, i.e. (t0, tf)
        Sol::Function   # Function providing the analytic solution of the IVP as a column vector
    end


    # ----- IVP1: Scalar exponential decay -----

       # Right-hand side of the ODE
       RHS1(t, y) = @. cos(π/4 + t) - (1e6)*(y - sin(π/4 +t))

       # Initial condition of the ODE
       IC1 = [0.707];
       
       # Integration timespan
       tspan1 = (0., 1. /10.);

        # Analytical solution of the ODE
        Sol1(t) = exp(-t);

        # Structure encoding the IVP
        IVP1 = InitialValueProblem(RHS1, IC1, tspan1, Sol1);


    # ----- IVP2: Harmonic oscillator with unitary amplitude and frequency -----

        # Right-hand side of the ODE
        ω = sqrt(2)
        ξ = 1/10
        RHS2(t, y::Vector) = [0 1; -ω^2 -ξ] * y;

        # Initial condition of the ODE
        IC2 = [0., 1.];
        
        # Integration timespan
        tspan2 = (0., 5.);
        
        # Analytical solution of the ODE
        Sol2(t) = [ sin(t/(2π)), cos(t/(2π)) ];

        # Structure encoding the IVP
        IVP2 = InitialValueProblem(RHS2, IC2, tspan2, Sol2);

end