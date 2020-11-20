using Test
using LevelSet
using StaticArrays
using LinearAlgebra
using Plots

h = 0.01
x = y = -2:h:2
grid = CartesianGrid(x,y)
m,n  = size(grid)
ϕ    = [1 - x^2 - y^2 for x in grid.xrange, y in grid.yrange]
𝐮    = [SVector(1.,0.) for _ in grid.xrange, y in grid.yrange]
∇ϕ   = upwind(ϕ,𝐮,grid)
∇ϕₑ  = [SVector(-2*x,-2*y) for x in grid.xrange, y in grid.yrange]
ee   = ∇ϕ - ∇ϕₑ
@test norm(ee[3:end-2,3:end-2],Inf) < 5*h

ϕ¹   = similar(ϕ)
Δt   = 0.1*h
nmax = 500
for n in 1:nmax
    evolve!(ϕ¹,ϕ,𝐮,grid,Δt)
    ϕ .= ϕ¹
end
fig = contour(x,y,transpose(ϕ),levels=[0],label="t=$(nmax*Δt)",aspect_ratio=:equal);
ϕₑ = (x,y) -> 1 - (x-0.5)^2 - y^2
contour!(fig,x,y,ϕₑ,levels=[0],label="t=$(nmax*Δt)")
