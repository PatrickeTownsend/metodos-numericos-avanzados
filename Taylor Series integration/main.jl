using Plots
using LinearAlgebra
using OffsetArrays
using LaTeXStrings
using ProgressMeter
using Dates
using AstroLib
include("functions.jl")
prog = ProgressUnknown("Procesando...")

μ = 3.986004418 * 1e5
t0 = 0
tf = 2*24*3600
h0 = 1
steps = round(Int,((tf - t0) / h0))
k_max = 20
Th = 2.2e-4
m0 = 24
Isp = 800
g0 = 9.80665
ṁ = -T/(Isp*g0) 
τ = 1e-15
x0 = [6678,0,0]
v0 = [0,sqrt(μ/norm(x0)),0]

#---Sun position--#
jd= juldate(DateTime(2000, 4, 20)) #year/month/day
vector = xyz(jd,2000)
r_sun=[vector[1]; vector[2]; vector[3]]

#-------SET ARRAYS----------#
X = [zeros(11)]
xn = OffsetArray(X[end], 0:10)
X = initializeX(xn, m0, x0, v0)
t = [t0]
iter = 0
step = [iter]
t_step = [h0]
thrust = [T]
#-------INITIALIZE----------#
while t[end] ≤ tf
    ProgressMeter.next!(prog)
    step = [step step[end]+1]
    xn = OffsetArray(X[:,end], 0:10)

    r_sc = [xn[0]; xn[1]; xn[2]]
    T = Th*Step_eclipse(r_sc,r_sun)
    thrust = [thrust T]
    global U = OffsetArray(zeros(11,k_max), 0:(11-1), 0:(k_max-1))
    global W = OffsetArray(zeros(14,k_max), 0:(14-1), 0:(k_max-1))
    U = initializeU(xn, U, μ, T, ṁ)
    W = initializeW(xn, U[:,0], W)
    U, W = Recursive(xn, U, W, k_max)
    uk1 = norm(U[:,k_max-2]./k_max)
    uk2 = norm(U[:,k_max-1]./k_max)
    h_next = 0.9*exp((1/(k_max-1))*log(τ/(uk1+h0*uk2)))
    t_step= [t_step h_next]
    tnext = t[end] + h_next
    X_next = xn + sum((U[:,k-1]/k)*(h_next)^k for k=1:k_max)
    X = [X X_next]
    t = [t tnext]
    h0 = h_next
    if iter ≥ 500000
        ProgressMeter.finish!(prog)
        break
    end
    sleep(0.1)
end
ProgressMeter.finish!(prog)
r0 = @. sqrt(X[1,1]^2 + X[2,1]^2 + X[3,1]^2)
rf = @. sqrt(X[1,end]^2 + X[2,end]^2 + X[3,end]^2)
println(size(thrust))
println(size(t))
plot(t[1,:],thrust[1,:])
plot(t' ./(24*3600),[X[1,:], X[2,:], X[3,:]], label=["X" "Y" "Z"], fontfamily="Computer modern", 
framestyle=:box, grid=true,tickfontsize=11,legendfontsize=11,titlefontsize=14,labelfontsize=12)
#yticks!([-8000,-6000,-4000,-2000, 0, 2000,4000,6000, 8000])
title!("Posición vs Tiempo")
xlabel!("Tiempo [Días]")
p1 = ylabel!("Posición [Km]")
#savefig(p1,"Taylor Series integration/plots/posicion.pdf")

plot(t' ./(24*3600),[X[4,:], X[5,:], X[6,:]], label=["Vx" "Vy" "Vz"], fontfamily="Computer modern", 
framestyle=:box, grid=true,tickfontsize=11,legendfontsize=11,titlefontsize=14,labelfontsize=12)
#yticks!([-8,-6,-4,-2, 0, 2,4,6, 8])
title!("Velocidad vs Tiempo")
xlabel!("Tiempo [Días]")
p2 = ylabel!("Velocidad [Km]")
#savefig(p2,"Taylor Series integration/plots/velocidad.pdf")

plot(X[1,:]./1e3,X[2,:]/.1e3, fontfamily="Computer modern", 
framestyle=:box, grid=true,tickfontsize=11,legendfontsize=11,titlefontsize=14,labelfontsize=12
, label="Orbita")
scatter!([0],[0],markersize=5, label="Planeta")
ylabel!(L"Y [km] $ x10^3$ ")
xlabel!(L"X [km] $ x10^3$ ")
p3 = title!("Orbita descrita")
#savefig(p3,"Taylor Series integration/plots/orbita.pdf")

p4 = plot(t'./(24*3600),X[7,:],title="Masa consumida", xlabel="Tiempo [Dias]", ylabel="Masa [Kg]",
fontfamily="Computer modern", 
framestyle=:box, grid=true,tickfontsize=11,legendfontsize=11,titlefontsize=14,labelfontsize=12,label=false)
#savefig(p4,"Taylor Series integration/plots/masa.pdf")

plot(step[2:end], t_step[2:end], title="Evolución del paso temporal", fontfamily="Computer modern", 
framestyle=:box, grid=true,tickfontsize=11,legendfontsize=11,titlefontsize=14,labelfontsize=12,legend=false)
xlabel!("Iteración")
xticks!([0,100,200,300,400,500,600,700,800])
yticks!([100,105,110,115,120,125,130,135,140,145])
p5= ylabel!(L"$Δt$ [s]")
#savefig(p5,"Taylor Series integration/plots/time_step.pdf")

display(p1)
display(p2)
display(p3)
display(p4)


