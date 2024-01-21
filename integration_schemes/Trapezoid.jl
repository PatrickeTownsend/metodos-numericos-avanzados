using MTH229

function TrapzSimple(f::Function, a::Float64, b::Float64)

    I = (b-a)*(f(b)+f(a))/2
    return I
end

function TrapzComp(f::Function, a::Float64, b::Float64, N::Int)
    h = (b-a)/N
    I = h*((f(a) + f(b))*0.5 + sum(f(a + k*h) for k=1:(N-1)))

    return I
end


f(x) = @. x^2

I = TrapzSimple(f,0,1)
I2 = TrapzComp(f, -1.0, 1.0, 1000)
I3 =  riemann(f,-1,1,1000, method="trapezoid")
display(I)
display(I2)
display(I3)
