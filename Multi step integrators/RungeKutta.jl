module RungeKutta

    #using ProgressBars

    export Explicit

# ----- Drivers for Explicit Runge-Kutta integrators -----

    # Driver for Explicit Runge-Kutta integrators
    function Explicit(IVP, order::Int, h)

        # Rename IVP variable for convenience
        t0 = IVP.tspan[1];
        tf = IVP.tspan[2];
        y0 = IVP.IC;

        # Numerical solution will be computed as we integrate step by step,
        # and stored in variables as follows:
        #  --  t = [t₀ t₁ t₂ t₃ ...]    where tᵢ are scalars representing time instants where the solution is provided
        #  --  y = [y₀ y₁ y₂ y₃ ...]    where yᵢ are column vectors representing the solution y(tᵢ) at times tᵢ
        #  --  f = [f₀ f₁ f₂ f₃ ...]    where fᵢ are column vectors representing the derivatives f(tᵢ, yᵢ)
        #
        # These variables must be first initialized, and the fill in as we integrate along our timespan.
        # We will do that appending the solutions tᵢ, yᵢ and fᵢ at each integration steps, although this is not particularly efficient.
        
        t = zeros(1);
        y = zeros(length(y0), 1);
        f = zeros(length(y0), 1);
        
        t[1]   =  t0;
        y[:,1] =  y0;
        f[:,1] =  IVP.RHS(t0, y0);
        
        # Now, take as many steps as need to cover the integration timespan.
        Steps = round(Int, (tf-t0) / h );

        # Create a ProgressBar object of prescribed length
        #pbar = ProgressBar(total=Steps)

        # Embed the 'Stepper' corresponding to the prescribed RK order
        if order == 1
            Stepper = RK1;
        
        elseif order == 2;
            Stepper = RK2;
        
        elseif order == 3;
            Stepper = RK3;
        
        elseif order == 4
            Stepper = RK4;
        
        end
        
        # Continue integration with the prescribed order 'Stepper'
        for i = 0:(Steps-1)

            # Take a step and evaluate the right-hand side
            ynext, fnext = Stepper(IVP.RHS, h, t[end], y[:,end], f[:,end]);
            
            # Update the solution at the new step
            t = [t (t0 + (i+1)*h)];
            y = [y ynext];
            f = [f fnext];

            # Update ProgressBar
            #update(pbar)

        end

        # When the integration is finished, return the solution as a tuple (t, y, f)
        return t, y, f

    end



# ----- Steppers for Explicit Runge-Kutta integrators -----

    # Stepper for Explicit Runge-Kutta of order 1: Forward Euler
    function RK1(RHS, h, tn, yn, fn)
        ynext = @. yn + h * fn
        fnext = RHS(tn, ynext)
        return ynext, fnext
    end


    # Stepper for Explicit Runge-Kutta of order 2: Heun's method
    function RK2(RHS, h, tn, yn, fn)
        a = [0 0; 1 0]
        b = [1/2 1/2]
        c = [0 1]
        
        k1 = fn
        k2 = RHS(tn + c[2]*h, yn .+ h.*a[2,1].*k1)
        
        ynext = @. yn + h * (b[1]*k1 + b[2]*k2)
        fnext = RHS(tn, ynext)
        
        return ynext, fnext
    end


    # Stepper for Explicit Runge-Kutta of order 3: Kutta's third-order method
    function RK3(RHS, h, tn, yn, fn)
        a = [0 0 0;
             1/2 0 0;
            -1 2 0]
        b = [1/6 2/3 1/6]
        c = [0 1/2 1]
        
        k1 = fn
        k2 = RHS( tn + c[2]*h, yn .+ h.*(a[2,1].*k1 ) )
        k3 = RHS( tn + c[3]*h, yn .+ h.*(a[3,1].*k1 .+ a[3,2].*k2) )
        
        ynext = @. yn + h * (b[1]*k1 + b[2]*k2 + b[3]*k3)
        fnext = RHS(tn, ynext)

        return ynext, fnext
    end


    # Stepper for Explicit Runge-Kutta of order 4: Classic fourth-order method
    function RK4(RHS, h, tn, yn, fn)
        a = [0 0 0 0;
             1/2 0 0 0;
             0 1/2 0 0;
             0 0 1 0]
        b = [1/6 1/3 1/3 1/6]
        c = [0 1/2 1/2 1]
        
        k1 = fn
        k2 = RHS( tn + c[2]*h, yn .+ h.*(a[2,1].*k1 ) )
        k3 = RHS( tn + c[3]*h, yn .+ h.*(a[3,1].*k1 .+ a[3,2].*k2) )
        k4 = RHS( tn + c[4]*h, yn .+ h.*(a[4,1].*k1 .+ a[4,2].*k2 .+ a[4,3].*k3) )
        
        ynext = @. yn + h * (b[1]*k1 + b[2]*k2 + b[3]*k3 + b[4]*k4)
        fnext = RHS(tn, ynext)

        return ynext, fnext
    end
    

end
