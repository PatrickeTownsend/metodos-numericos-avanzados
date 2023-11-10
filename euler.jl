using Plots

F(t,y)=-y

sol(t) = @. exp(-t)

function Euler(tn,yn,h)
    y_next = yn + h*F(tn,yn)
    return y_next
end
function BDF1(tn1,yn1,h)
    y_next = yn1 + h*F(tn1,yn1)
    return y_next
end
function BF2(tn,yn,h)
    ynext = -(4/3)yn + (1/3)yn_1 + (2/3)h*F(tn,yn)
    return ynext
end

y0 = 1.
t0 = 0.0
tf = 2
h = 0.01
ysol = y0
tsol = t0;
#---Explicit Euler---#
ycurrent = y0
for i=0:round(Int,(tf-t0)/h)-1
    ynext = Euler(t0+i*h,ycurrent,h)
    ysol = [ysol ynext]
    tsol = [tsol t0+(i+1)*h]
    ycurrent = ynext
end
#---BDF1----#
ycurrent = y0
for i=0:round(Int,(tf-t0)/h)-1
    ynext = BDF1(t0+(i+1)*h,ynext,h)
    ysol = [ysol ynext]
    tsol = [tsol t0+(i+1)*h]
    ycurrent = ynext
end
yexact = sol(tsol)
plot(tsol',abs.((ysol'-yexact')./yexact'))