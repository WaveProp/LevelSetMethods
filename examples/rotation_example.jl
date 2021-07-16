using Test
using LevelSetMethods
using WavePropBase
using LinearAlgebra
using Plots

nx,ny = 201,201
x     = LinRange(-2,2,nx)
y     = LinRange(-2,2,ny)
M  = UniformCartesianMesh(x,y)
bc    = PeriodicBC(3)
ϕ₀    = LevelSet(M,bc) do (x,y)
    sqrt((x+0.4)^2 + y^2) - 0.35
end
𝐮     = NodeField(M) do (x,y)
    2π*SVector(-y,x)
end
term1  = AdvectionTerm(velocity=𝐮,scheme=Upwind())
terms  = (term1,)
fig = plot(ϕ₀)

# Forward  euler
integrator = ForwardEuler()
eq = LevelSetEquation(;terms,integrator,state=deepcopy(ϕ₀),t=0,cfl=0.2)
@time integrate!(eq,1)
plot!(fig,eq,lc=:blue)

# Low storage RK integrator
integrator = RKLM2()
eq = LevelSetEquation(;terms,integrator,state=deepcopy(ϕ₀),t=0)
@time integrate!(eq,1)
plot!(fig,eq,lc=:red,m=:x)

# anim = @animate for n ∈ 0:100
#     tf = 0.02*n
#     integrate!(eq,tf)
#     plot(eq)
# end
# gif(anim, "test.gif")
