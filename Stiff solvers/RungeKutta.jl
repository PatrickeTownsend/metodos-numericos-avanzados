module RungeKutta
   export RK4

   function RK4(IVP, h)
    t0 = IVP.tspan[1]
    tf = IVP.tspan[2]
    y0 = IVP.IC

    Steps = round(Int, (tf-t0)/h)

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
    end
    return t, y
end


end