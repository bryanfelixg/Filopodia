#=
Looking at 
du/dt = D1(A_yy u + u A_xx) + au^2/v + ubar - au
dv/dt = D2(A_yy v + v A_xx) + au^2   - βv
=#
using LinearAlgebra
using DifferentialEquations
using BenchmarkTools
using Sundials
using Plots
# Some constants

p  = (1.0,1.0,1.0,10.0,.01,100.0) # a,α,ubar,β,D1,D2
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
r0 = zeros(100,100,2)
r0[:,:,1] .= uss.+0.1.*rand.()
r0[:,:,2] .= vss #.+1.0.*rand.()

## Problem
prob = ODEProblem(basic!,r0,(0.0,500.0),p)
sol = solve(prob,CVODE_BDF(linear_solver=:GMRES),
progress=true,
save_everystep=false,
saveat = [0.0, 100.0, 200.0, 300.0, 400.0, 500.0]) 

##
X = reshape([i for i in 1:100 for j in 1:100],N,N)
Y = reshape([j for i in 1:100 for j in 1:100],N,N)

p1 = surface(X,Y,sol.u[3][:,:,1],title = "[A]",camera=(0,90))
p2 = surface(X,Y,sol.u[3][:,:,2],title = "[B]",camera=(0,90))
plot(p1,p2,layout=grid(2,1)) 

##
animate(sol)