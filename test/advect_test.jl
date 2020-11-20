using Test
using LevelSetMethods
using LinearAlgebra

hx,hy = 0.01, 0.02
x = collect(-1:hx:1)
y = collect(-2:hy:2)
grid = CartesianGrid(x,y)
m,n  = size(grid)
ϕ    = LevelSet(grid) do (x,y)
    1 - x^2 - y^2
end    
𝐮    = MeshField(x->SVector(x[1],x[2]),grid)
adv = AdvectionTerm(velocity=𝐮)

buffer = adv(ϕ) 
ref  = [-2*x^2-2*y^2 for (x,y) in grid]
ee   = values(buffer) - ref
@test norm(ee[3:end-2,3:end-2],Inf) < 5*max(hx,hy)

buffer = zero(ϕ)
adv(buffer,ϕ) # shoud be allocation-free
ee   = values(buffer) - ref
@test norm(ee[3:end-2,3:end-2],Inf) < 5*max(hx,hy)

