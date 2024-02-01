using Plots
using LinearAlgebra
using OffsetArrays
include("functions.jl")

μ = 3.986004418 * 1e5
t0 = 0
tf = 3.821*24*3600
h = 10
steps = round(Int,((tf - t0) / h))
k_max = 20
T = 2e-4
m0 = 100
Isp = 3000
g0 = 9.80665
ṁ = -T/(Isp*g0) 

x0 = [6678,0,0]
v0 = [0,sqrt(μ/norm(x0)),0]

#-------SET ARRAYS----------#
X = OffsetArray(zeros(11,steps), 0:(11-1), 0:(steps-1))
X = initializeX(X, m0, x0, v0)
t = OffsetArray(zeros(steps), 0:(steps-1))
#-------INITIALIZE----------#
for i = 1:steps-1
    xn = OffsetArray(X[:,i-1], 0:10)
    t[i] = t[i-1] + h
    global U = OffsetArray(zeros(11,k_max), 0:(11-1), 0:(k_max-1))
    global W = OffsetArray(zeros(14,k_max), 0:(14-1), 0:(k_max-1))
    U = initializeU(xn, U, μ, T, ṁ)
    W = initializeW(xn, U[:,0], W)
    U, W = Recursive(X, U, W, k_max, i-1)
    X[:,i] = X[:,i-1] + sum((U[:,k-1]/k)*(t[i]-t[i-1])^k for k=1:k_max)
end

plot(t,[X[0,:], X[1,:], X[2,:]], label=["X" "Y" "Z"])
yticks!([-8000,-6000,-4000,-2000, 0, 2000,4000,6000, 8000])

plot(t,[X[3,:], X[4,:], X[5,:]], label=["Vx" "Vy" "Vz"])
yticks!([-8,-6,-4,-2, 0, 2,4,6, 8])
println((X[0,0:10]))



