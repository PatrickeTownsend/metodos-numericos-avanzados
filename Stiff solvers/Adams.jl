module Adams
    export Bashforth
    using LinearAlgebra
    include("LUdescomp.jl")


    function Bashforth(IVP, k::Int, h)
        t0 = IVP.tspan[1]
        tf = IVP.tspan[2]
        y0 = IVP.IC
        function_eval = 0

        Steps = round(Int, (tf-t0)/h)
        println(Steps)

        Current_k = 1
        Stepper = AB1

        t = zeros(Steps)
        y = zeros(length(y0), Steps)
        f = zeros(length(y0), Steps)
        y[:,1] += y0
        t[1] += t0
        f[:,1] += IVP.RHS(t0, y0)
        function_eval += 1


        for n = 1:(k-1)

            t[n+1] += t[n] + h
            y[:,n+1] += Stepper(h, y[:,n], f[:,1:n])
            f[:,n+1] += IVP.RHS(t[n+1], y[:,n+1])
            function_eval += 1
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
            y[:,n+1] += Stepper(h, y[:,n], f[:,n-(k-1):n])
            f[:,n+1] += IVP.RHS(t[n], y[:,n])
            function_eval += 1
            t[n+1] += t[n] + h
        end

        return t, y, f, function_eval
    end

    function Moulton(IVP,k::Int, h, prob)
        t0 = IVP.tspan[1]
        tf = IVP.tspan[2]
        y0 = IVP.IC
        tol = 1e-8
        max_iter = 5000
        function_eval = 0
        iter_global = 0

        Steps = round(Int, (tf-t0)/h)
        println(Steps)

        Current_k = 1
        Stepper = AB1

        t = zeros(Steps)
        y = zeros(length(y0), Steps)
        f = zeros(length(y0), Steps)
        y[:,1] += y0
        t[1] += t0
        f[:,1] += IVP.RHS(t0, y0)
        function_eval += 1
        Jac_eval = 0


        for n = 1:(k-1)

            t[n+1] += t[n] + h
            y[:,n+1] += Stepper(h, y[:,n], f[:,1:n])
            f[:,n+1] += IVP.RHS(t[n+1], y[:,n+1])
            function_eval += 1
            Current_k += 1

            if  Current_k == 2
                Stepper = AB2
            elseif Current_k == 3
                Stepper = AB3
            elseif Current_k == 4
                Stepper = AB4
            end
        end

        if prob == 2
            for n = k:(Steps-1)
                t[n+1] += t[n] + h
                δ = 1
                y[:,n+1] = y[:,n]
                m = size(y)[1]

                for i= 1:max_iter
                    J = Matrix{Float64}(I,m,m) .- h .*(9/24).*IVP.Jacobian(y[:,n+1])
                    Jac_eval+=1
                    iter_global += 1
                    f = -y[:,n] - (h/24)*(9*IVP.RHS(t[n+1],y[:,n+1]) + 19*IVP.RHS(t[n],y[:,n]) -5*IVP.RHS(t[n-1],y[:,n-1])
                    + IVP.RHS(t[n-2],y[:,n-2])) + y[:,n+1]
                    δ = - J \ f
                    function_eval += 1
                    if norm(δ) ≤ tol
                        println("Newton converged at $i iterations",y[:,n+1])
                        break
                    end
                    y[:,n+1] += δ
                end
                
                
            end
        else
            for n = k:(Steps-1)
                t[n+1] += t[n] + h
                ynext, fnext = AM3(h, y[:,n], f[:,n-2:n], t[n+1], IVP, tol)
                y[:,n+1] += ynext
                f[:,n+1] += fnext
            end
        end

        println(iter_global)


        return t, y, function_eval, Jac_eval

    end

    function AM3(h, yn, f, t, IVP, tol)
        δ = 1
        ynext  =  yn
        fnext = f[:,end]
        m = size(yn)[1]
        iter = 0
        while norm(δ)>tol
            J = Matrix{Float64}(I, m,m) .- h .*(9/24).*IVP.Jacobian(ynext)
            L,U = LUfactor(J)
            Fₙ = -yn - (h/24)*(9*fnext + 19*f[:,end] -5*f[:,end-1]
                    + f[:,end-2]) + ynext
            δ = LU_linear_solve(L,U,Fₙ)
            ynext += δ
            fnext += IVP.RHS(t,ynext)
            iter += 1
        end
        println("Newton converged at $iter iterations")
        return ynext, fnext
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