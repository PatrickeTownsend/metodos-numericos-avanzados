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
xₗ  = ChebyshevNode(n) #LinRange(x_init, x_end, n)
yᵢ = func(xᵢ)
yₗ  = func(xₗ)
Pᵢ = zeros(n)
Iᵢ = zeros(n)
ϵ = zeros(n)
for i in eachindex(xₗ)
    Pᵢ[i] += Lagrange(xₗ[i] ,xᵢ ,yᵢ)[1]
    ϵ[i] += epsilon(xₗ[i], xᵢ)
end
Eₗ = @.abs(yₗ - Pᵢ) 
#Plot_func(xₗ, yₗ, Pᵢ) #Grafica la funciones exacta e interpolada

scatter(xₗ , [log.(Eₗ), log.(abs.(ϵ))], title="Error de interpolacion", xlabel="xₗ", ylabel="log(Error)", label=["|f(x) -F(x)|" "ϵ(x)"],grid=true)





