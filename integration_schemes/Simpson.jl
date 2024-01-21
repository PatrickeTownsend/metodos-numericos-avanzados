using MTH229

function Simpson_13_Simp(f::Function, a::Float64, b::Float64)
    h = (b-a)/2
    m = (a+b)/2
    I = (h/3)*(f(a) + 4*f(m) + f(b))
    return I
end

function Simpson_13_Comp(f::Function, a::Float64, b::Float64, N::Int)
    h = (b-a)/N
    x = zeros(N)
    for i=1:N
        x[i] = f(a + i*h)
    end
    upper1 = round(Int,(0.5*N-1))
    upper2 = round(Int,(N*0.5))
    I = (h/3)*(f(a) + 2*sum(x[2*k] for k=1:upper1) + 4*sum(x[2*k-1] for k=1:upper2) + x[N])
    return I
end


function Simpson_38_Simp(f::Function, a::Float64, b::Float64)
    h = (b - a)/3
    I = (3/8)*h*(f(a) + 3*f((2*a+b)/3) + 3*f((a+2*b)/3) + f(b))
    return I
end

function Simpson_38_Comp(f::Function, a::Float64, b::Float64, N::Int)
    h = (b-a)/N
    x = zeros(N)
    for i=1:N
        x[i] = f(a + i*h)
    end
    upper1 = round(Int,(N/3 -1))
    upper2 = round(Int,(N/3 -2))


    I = (3/8)*h*(x[1] + 3*sum(x[3*i+1] for i=0:upper1) + 3*sum(x[3*i+2] for i=0:upper1) + 
    2*sum(x[3*i+3] for i=0:upper2) + x[N])

    return I
end

f(x) = @. x^2 + x^3
I = Simpson_13_Simp(f,0.0, 1.0)
I2 = Simpson_13_Comp(f, 0.0, 1.0, 1000)
I3 = Simpson_38_Simp(f, 0.0, 1.0)
I4 = Simpson_38_Comp(f, 0.0, 1.0, 3000)
I5 = riemann(f,0,8,1000, method="simpsons")
display(I)
display(I2)
display(I3)
display(I4)