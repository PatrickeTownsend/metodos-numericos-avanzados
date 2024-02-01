module RungeKutta
   export RK4
   export SDIRK4
   export SDIRK4_v2
   using LinearAlgebra

   function RK4(IVP, h)
      t0 = IVP.tspan[1]
      tf = IVP.tspan[2]
      y0 = IVP.IC
      function_eval = 0

      Steps = round(Int, (tf-t0)/h)
      println(Steps)

      t = zeros(Steps)
      y = zeros(length(y0), Steps)
      f = zeros(length(y0), Steps)
      y[:,1] += y0
      t[1] += t0

      for n=1:(Steps-1)
          k1 = IVP.RHS(t[n], y[:,n])
          k2 = IVP.RHS(t[n] + h/2, y[:,n] + k1 * h/2)
          k3 = IVP.RHS(t[n] + h/2, y[:,n] + k2 * h/2)
          k4 = IVP.RHS(t[n] + h, y[:,n] + k3 * h)
          y[:,n+1] += y[:,n] + (h/6)*(k1 + 2*k2 +2*k3 + k4)
          t[n+1] += t[n] + h
          function_eval += 4
       end
       return t, y, function_eval
    end


    function SDIRK4(IVP, h)
        A, b, c = Butcher5s()
        t0 = IVP.tspan[1]
        tf = IVP.tspan[2]
        y0 = IVP.IC

        Steps = round(Int, (tf-t0)/h)
        println(Steps)

        t = zeros(Steps)
        y = zeros(length(y0), Steps)
        f = zeros(length(y0), Steps)
        
        y[:,1] += y0
        t[1] += t0

        d = length(y0)
        s = length(b)

        #----Newton-Raphson set-up----#

        max_tol = 1e-16
        max_iter = 5000
        κ = 1e-2
        Iden = Matrix{Float64}(I,d,d)

        Mat = kron(A,Iden)  # Matriz sd x sd
    
        di = A \ b
        println(size(di))
        function_eval = 0
        Jac_eval = 0
        iter_global = 0

        for n=1:(Steps-1)
            zk = zeros(s*d)

            tn = t[n]
            yn = y[:,n]
            J = IVP.Jacobian(yn)
            Jac_eval +=1
            function_eval += 1
            JJ = Matrix{Float64}(I,s*d, s*d) - h.*kron(A,J)
            iter = 0

            delta_z = ones(2)
            ηₖ = 1

            while ηₖ*delta_z[2] ≥ κ*max_tol
                iter += 1
                iter_global += 1
                delta_z[1] = delta_z[2]
                Fzᵏ= [h*IVP.RHS(tn+c[1]*h, yn + zk[d*1-d+1:d*1]);
                      h*IVP.RHS(tn+c[2]*h, yn + zk[d*2-d+1:d*2]);
                      h*IVP.RHS(tn+c[3]*h, yn + zk[d*3-d+1:d*3]);
                      h*IVP.RHS(tn+c[4]*h, yn + zk[d*4-d+1:d*4]);
                      h*IVP.RHS(tn+c[5]*h, yn + zk[d*5-d+1:d*5])]
                function_eval += 5
                F = -zk + h.*Mat*Fzᵏ
                Δzᵏ= -JJ \ F
                delta_z[2] = norm(Δzᵏ)
                zk +=  Δzᵏ
                θₖ = delta_z[2]/delta_z[1]
                ηₖ = θₖ/(1-θₖ)
                #println(delta_z[2])
                if iter > 10000
                    break
                end
                
            end
            zᵢ=[zk[d*1-d+1:d*1];
                zk[d*2-d+1:d*2];
                zk[d*3-d+1:d*3];
                zk[d*4-d+1:d*4];
                zk[d*5-d+1:d*5]]

            t[n+1]= tn + h
            phi = zeros(d)
            for i= 1:s
                phi = phi + di[i]*zk[d*i-d+1:d*i]
            end
            y[:,n+1] =  yn + phi
            println("Newton converged with $iter iterations",y[:,n+1])

        end
        println(iter_global)
        return t, y, function_eval, Jac_eval

    end

    function Butcher5s()
        A = [1/4  0  0 0 0;
             1/2 1/4 0 0 0;
             17/50 -1/25 1/4 0 0;
             371/1630 -137/2720 15/544 1/4 0;
             25/24 -49/48 125/16 -85/12 1/4]
        
        b = [25/24, -49/48, 125/16, -85/12, 1/4]

        c = [1/4 3/4 11/20 1/2 1]
        return A, b, c
    end

end

