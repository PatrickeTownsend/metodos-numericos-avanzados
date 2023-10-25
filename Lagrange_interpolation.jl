using PlotlyJS
using Plots
include("interp_functions.jl")

#---Funcion---#
func(x) = @.  sin(x) #1/(1+25*x^2)
#---Intervalo---#
x_init = -10.0
x_end = 10.0

#---Malla de puntos---#
N = 20 # Numero de nodos interpolantes
n = 1000 # Numero de puntos a graficar

#---Arrays---#
xᵢ = LinRange(x_init, x_end, N)
xₗ  =  LinRange(x_init, x_end, n) #ChebyshevNode(x_init,x_end,n) #
yᵢ = func(xᵢ)
yₗ  = func(xₗ)
Pᵢ = zeros(n)
Iᵢ = zeros(n)
ϵ = zeros(n)
for i in eachindex(xₗ)
    Pᵢ[i] += Lagrange(xₗ[i] ,xᵢ ,yᵢ)[1]
    ϵ[i] += epsilon_eq(xₗ[i], xᵢ)
end
Eₗ = @.abs(yₗ - Pᵢ) 
Plot_func(xₗ, yₗ, Pᵢ) #Grafica la funciones exacta e interpolada

Plots.scatter(xₗ , [log.(Eₗ), log.(abs.(ϵ))], title="Error de interpolacion chevyshev", xlabel="xₗ", ylabel="log(Error)", label=["|f(x) -F(x)|" "ϵ(x)"],grid=true)

trace_error = PlotlyJS.scatter(x=xₗ, y=Eₗ, mode="markers", name="|f(x)-F(x)|")
trace_epsilon = PlotlyJS.scatter(x=xₗ, y=ϵ, mode="markers", name="ϵ(x)")
layout = PlotlyJS.Layout(title="error de interpolacion chebyshev",
                xaxis_title="x",
                yaxis_title="error",
                yaxis_type="log")  

PlotlyJS.plot([trace_error, trace_epsilon], layout)
