using Plots
include("interp_functions.jl")

#---Funcion---#
func(x) = @. 1/(1+25*x^2)
#---Intervalo---#
x_init = -10.0
x_end = 10.0

#---Malla de puntos---#
N = 20 # Numero de nodos interpolantes
n = 1000 # Numero de puntos a graficar

#---Arrays---#
xᵢ = LinRange(x_init, x_end, N)
xₗ  = LinRange(x_init, x_end, n)
yᵢ = func(xᵢ)
Pᵢ = zeros(N)


