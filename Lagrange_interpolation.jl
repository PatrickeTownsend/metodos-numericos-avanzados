using PlotlyJS
using Plots
include("interp_functions.jl")

#---Funcion---#
func(x) = @.  1/(1+25*x^2)
#func(x) = @. exp(-(5*x)^2)
#---Intervalo---#
x_init = -1.0
x_end = 1.0

#---Malla de puntos---#
N = 20 # Numero de nodos interpolantes
n = 1000 # Numero de puntos a graficar
dist = "B"
funcion = "1"
#---Arrays---#
a = LinRange(x_init,x_end,N)
#xᵢ = LinRange(x_init, x_end, N)
xᵢ = distribucionB(x_init,x_end,N)
#xᵢ = distribucionC(x_init,x_end,N)
xₗ  = LinRange(x_init, x_end, n) 
yᵢ = func(xᵢ)
yₗ = func(xₗ)
Pᵢ = zeros(n)
Iᵢ = zeros(n)
ϵ = zeros(n)

for i in eachindex(xₗ)
    Pᵢ[i] += Lagrange(xₗ[i] ,xᵢ ,yᵢ)
    ϵ[i] += epsilon(xₗ[i], xᵢ)
   
end
#---Bruja---#
Eₗ = @.abs(yₗ - Pᵢ) 
Plot_func(xₗ, yₗ, Pᵢ) #Grafica la funciones exacta e interpolada

#---Funcion 1 equispaced---#
trace_error = PlotlyJS.scatter(x=xₗ, y=Eₗ, mode="lines", name="|f(x)-F(x)|")
trace_epsilon = PlotlyJS.scatter(x=xₗ, y=abs.(ϵ), mode="lines", name="ϵ(x)")
layout = PlotlyJS.Layout(title="error de interpolacion |ϵ(x)| dist $dist funcion $funcion",
                xaxis_title="x",
                yaxis_title="error |ϵ(x)|",
                yaxis_type="log")  

layout2 = PlotlyJS.Layout(title="error de interpolacion f₁(x) dist $dist",
                xaxis_title="x",
                yaxis_title="error",
                yaxis_type="log")  
interp_error = PlotlyJS.plot([trace_error, trace_epsilon], layout2)
error_plot = PlotlyJS.plot([trace_epsilon], layout)
display(interp_error)
display(error_plot)

PlotlyJS.savefig(error_plot, "epsilon_$dist $funcion.png")
PlotlyJS.savefig(interp_error, "error_tot_$dist $funcion.png")

