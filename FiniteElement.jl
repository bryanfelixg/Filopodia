using JuAFEM, SparseArrays

## MESH
grid = generate_grid(Quadrilateral,(100,100))


## Test Functions
dim = 2
ip = Lagrange{dim,RefCube,1}()
qr = QuadratureRule{dim, RefCube}(2)
cellvalues = CellScalarValues(qr,ip)


## Degrees of Freedom
dh = DofHandler(grid)
push!(dh, :u, 1)
close!(dh)

## Sparsity
K = create_sparsity_pattern(dh)
M = create_sparsity_pattern(dh)

f = zeros(ndofs(dh)) # Preallocation

##Boundary Conditions

max_temp = 100
Δt = 1
T = 200
ch = ConstraintHandler(dh)

∂Ω1 = union(getfaceset.((grid, ), ["left","right"])...)
dbc = Dirichlet(:u, ∂Ω1, (x,t)->0)
add!(ch,dbc)

∂Ω2 = union(getfaceset.((grid, ), ["top","bottom"])...)
dbc = Dirichlet(:u, ∂Ω2, (x,t)->t*(max_temp/T))
add!(ch,dbc)
close!(ch)
update!(ch,0.0)

include("doassemble.jl")

## Solution
K, f = doassemble_K!(K,f,cellvalues,dh)
M = doassemble_M!(M,cellvalues,dh)
A = (Δt .* K) + M

## RHS

rhsdata = get_rhs_data(ch,A)

uₙ = 5.0.*ones(length(f))
apply!(A,ch)

pvd = paraview_collection("results/transient-heat.pvd")

##Solve


for t in 0:Δt:T
    update!(ch,t)

    b = Δt .* f.+ M * uₙ
    apply_rhs!(rhsdata, b, ch)
    
    u = A \ b

    vtk_grid("results/transient-heat-$t",dh) do vtk
        vtk_point_data(vtk,dh,u)
        vtk_save(vtk)
        pvd[t]=vtk
    end

    uₙ .= u
end
