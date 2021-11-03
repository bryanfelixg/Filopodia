#=
Looking at 
du/dt = D1(A_yy u + u A_xx) + au^2/v + ubar - au
dv/dt = D2(A_yy v + v A_xx) + au^2   - βv
=#
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using BenchmarkTools
using Sundials
using Plots
## Some constants
tspan = (0.0,100.0)
p  = (1.0,1.0,1.0,10.0,10,1.0) # a,α,ubar,β,D1,D2
N  = 100
Axx = Array(Tridiagonal(
    [1.0 for i in 1:N-1],
    [-2.0 for i in 1:N],
    [1.0 for i in 1:N-1]
    ))
Ayy = copy(Axx)
Axx[2,1] = 2.0 
Axx[end-1,end] = 2.0
Ayy[1,2] = 2.0
Ayy[end,end-1] = 2.0

##
#Some allocation
Ayyu = zeros(N,N)
uAxx = zeros(N,N)
Du = zeros(N,N)
Ayyv = copy(Ayyu)
vAxx = copy(uAxx)
Dv = copy(Du)
function basic!(dr,r,p,t)
    a,α,ubar,β,D1,D2 = p
    u = @view r[:,:,1]
    v = @view r[:,:,2]
    du = @view dr[:,:,1] 
    dv = @view dr[:,:,2]
    mul!(Ayyu,Ayy,u)
    mul!(uAxx,u,Axx)
    mul!(Ayyv,Ayy,v)
    mul!(vAxx,v,Axx)
    @. Du = D1*(Ayyu + uAxx)
    @. Dv = D2*(Ayyv + vAxx)
    @. du = Du + a*u*u/v + ubar - α*u
    @. dv = Dv + a*u*u   - β*v
end

## Initial Conditions
a,α,ubar,β,D1,D2 = p
uss = (ubar +β)/α
vss = (a/β)*uss^2
r0 = zeros(N,N,2)
r0[:,:,1] .= 0 #uss.+0.1*rand.()
r0[:,:,2] .= 0 #vss 

include("dos.jl")
## Problem
prob = ODEProblem(basic!,r0,tspan,p)
sol = solve(prob,CVODE_BDF(linear_solver=:GMRES),
progress=true,
save_everystep=false,
saveat = range(tspan[1],tspan[2], length=100)) 

##
clims = (0,50);
for i=1:size(sol.t)[1]
    gr()
    data1 = @view sol.u[i][:,:,1]
    data2 = @view sol.u[i][:,:,2]
    p1=heatmap(1:size(data1,1),
        1:size(data1,2), data1,
        c=cgrad([:blue, :white,:red, :yellow]),
        xlabel="x values", ylabel="y values",
        clims = clims,
        title="[A]")
    p2=heatmap(1:size(data2,1),
        1:size(data2,2), data2,
        c=cgrad([:blue, :white,:red, :yellow]),
        xlabel="x values", ylabel="y values",
        clims = clims,
        title="[A]")
    plot(p1,p2,layout=grid(2,1)) |>display
end


##
# X = reshape([i for i in 1:100 for j in 1:100],N,N)
# Y = reshape([j for i in 1:100 for j in 1:100],N,N)
# p1 = surface(X,Y,sol.u[3][:,:,1],title = "[A]",camera=(0,90))
# p2 = surface(X,Y,sol.u[3][:,:,2],title = "[B]",camera=(0,90))
# plot(p1,p2,layout=grid(2,1)) 

####animate(sol)
