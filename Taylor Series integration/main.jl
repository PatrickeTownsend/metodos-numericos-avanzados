using Plots
using LinearAlgebra
using OffsetArrays
include("functions.jl")

μ = 3.986004418 * 1e5
t0 = 0
tf = 3.821*24*3600
h = 10
steps = round(Int,((tf - t0) / h))
k = 20
T = 2e-4
m0 = 100
Isp = 3000
g0 = 9.80665
ṁ = -T/(Isp*g0) 

x0 = [6678,0,0]
v0 = [0,sqrt(μ/norm(x0)),0]

#-------SET ARRAYS----------#
X = OffsetArray(zeros(11,k), 0:(11-1), 0:(Steps-1))
U = OffsetArray(zeros(11,k), 0:(11-1), 0:(k-1))
W = OffsetArray(zeros(14,k), 0:(14-1), 0:(k-1))

#-------INITIALIZE----------#
X = initializeX(X, m0, x0, v0)
U = initializeU(X, U, μ, T, ṁ)
W = initializeW(X, U, W)

