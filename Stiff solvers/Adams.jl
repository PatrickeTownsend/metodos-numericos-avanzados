module Adams
    export Bashforth


    function Bashforth(IVP, k::Int, h)
        t0 = IVP.tspan[1]
        tf = IVP.tspan[2]
        y0 = IVP.IC

        Steps = round(Int, (tf-t0)/h)

        Current_k = 1
        Stepper = AB1

        t = zeros(Steps)
        y = zeros(length(y0), Steps)
        f = zeros(length(y0), Steps)
        y[:,1] += y0
        t[1] += t0
        f[:,1] += IVP.RHS(t0, y0)


        for n = 1:(k-1)

            t[n+1] += t[n] + h
            y[:,n+1] += Stepper(h, y[:,n], f[:,1:n])
            f[:,n+1] += IVP.RHS(t[n+1], y[:,n+1])
            Current_k += 1

            if  Current_k == 2
                Stepper = AB2
            elseif Current_k == 3
                Stepper = AB3
            elseif Current_k == 4
                Stepper = AB4
            end
        end


        for n = (k):(Steps-1)
            y[:,n+1] += Stepper(h, y[:,n], f[:,n-3:n])
            f[:,n+1] += IVP.RHS(t[n+1], y[:,n+1])
            t[n+1] += t[n] + h
        end
        
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