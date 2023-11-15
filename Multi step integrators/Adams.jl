module Adams

    export Bashforth

# ----- Drivers for Adams integrators -----

    # Driver for Adams-Bashforth integrators
    function Bashforth(IVP, k::Int, h)

        println("Computing numerical solution...")
        
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

        # Start-up procedure
        Current_k = 1;
        Stepper = AB1;
        
        for i = 0:(k-2)

            # 'i' defines, during the start-up, how may points behind the current time
            # must be passed to the Stepper at the current step

            # Take a step and evaluate the right-hand side
            ynext = Stepper(h, y[:,end], f[:, end-i:end]);
            fnext = IVP.RHS(t[end], y[:,end])

            # Update the solution at the new step
            t = [t (t0 + (i+1)*h)];
            y = [y ynext];
            f = [f fnext];

            # Use appropriate stepper, so at each step we increase
            # the order of the method as past point become available.
            Current_k += 1;
            
            if Current_k == 2
                Stepper = AB2;
            elseif Current_k == 3
                Stepper = AB3;
            elseif Current_k == 4
                Stepper = AB4;
            elseif Current_k == 5
                Stepper = AB5;
            end
        
        end
        
        # Continue integration with the prescribed order 'Stepper'
        for i = (k-1):(Steps-1)

            ynext = Stepper(h, y[:,end], f[:, end-k+1:end]);
            fnext = IVP.RHS(t[end], y[:,end])

            t = [t (t0 + (i+1)*h)];
            y = [y ynext];
            f = [f fnext];

        end

        # When the integration is finished, return the solution as a tuple (t, y, f)
        return t, y, f

    end



# ----- Steppers for Adams-Bashforth integrators -----

    # Stepper for Adams-Bashforth of order 1
    function AB1(h, yn, fn)
        ynext = @. yn + h * fn
        return ynext
    end


    # Stepper for Adams-Bashforth of order 2
    function AB2(h, yn, fprev)
        # Here the input 'fprev' contain the right-hand side evaluation not only in the current time, tn,
        # but also in the preceding instant of time, tn-1.
        #  --  fpref = [f(tn-1) f(tn)]
        ynext = @. yn + h * 0.5( 3. * fprev[:,end] - fprev[:,end-1] );
        return ynext
    end


    # Stepper for Adams-Bashforth of order 3
    function AB3(h, yn, fprev)
        # Here the input 'fprev' contain the right-hand side evaluation not only in the current time, tn,
        # but also in the preceding instants of time.
        #  --  fpref = [f(tn-2) f(tn-1) f(tn)]
        ynext = @. yn + h * ( 23. * fprev[:,end] 
                            - 16. * fprev[:,end-1] 
                            +  5. * fprev[:,end-2] ) / 12.;
        return ynext
    end


    # Stepper for Adams-Bashforth of order 4
    function AB4(h, yn, fprev)
        # Here the input 'fprev' contain the right-hand side evaluation not only in the current time, tn,
        # but also in the preceding instants of time.
        #  --  fpref = [f(tn-3) f(tn-2) f(tn-1) f(tn)]
        ynext = @. yn + h * ( 55. * fprev[:,end] 
                            - 59. * fprev[:,end-1] 
                            + 37. * fprev[:,end-2] 
                            -  9. * fprev[:,end-3] ) / 24.;
        return ynext
    end


    # Stepper for Adams-Bashforth of order 5
    function AB5(h, yn, fprev)
        # Here the input 'fprev' contain the right-hand side evaluation not only in the current time, tn,
        # but also in the preceding instants of time.
        #  --  fpref = [f(tn-4) f(tn-3) f(tn-2) f(tn-1) f(tn)]
        ynext = @. yn + h * ( 1901. * fprev[:,end] 
                            - 2774. * fprev[:,end-1] 
                            + 2616. * fprev[:,end-2] 
                            - 1274. * fprev[:,end-3] 
                            +  251. * fprev[:,end-4] ) / 720.;
        return ynext
    end
    

end
