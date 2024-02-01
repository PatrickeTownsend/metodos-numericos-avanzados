#-----Patrick Townsend-----#
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

#----constantes----#
ω = sqrt(2)
ξ = 1/10
h = 0.01
t0 = 0
tf = 5.0
y0 = [0. ,1.]

#----Set-up----#
RHS(t,y) = [0 1;-ω^2 -ξ]* y
Steps = round(Int,(tf-t0)/h)
t = zeros(Steps)
y = zeros(length(y0), Steps)
f = zeros(length(y0), Steps)
y[:,1] += y0
t[1] += t0

for n=1:(Steps-1)
    ynext, fnext = RK3(RHS, h, t[n], y[:,n], f[:,n])
    y[:,n+1] += ynext
    f[:,n+1] += fnext
    t[n+1] += t0 + n*h
end

display(y[:,end])
display(t[end])


