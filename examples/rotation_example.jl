using Test
using LevelSetMethods
using LinearAlgebra
using StaticArrays
using Plots

import WavePropBase as WPB

nx,ny = 200,200
rec   = WPB.HyperRectangle((-2,-2),(2.,2.))
M     = WPB.UniformCartesianMesh(rec,(nx,ny))
bc    = PeriodicBC(3)
ϕ     = DiscreteLevelSet(M,0) do (x,y)
    sqrt((x+0.4)^2 + y^2) - 0.35
end
𝐮     = CartesianGridFunction(M) do (x,y)
    2π*SVector(-y,x)
end
term1  = AdvectionTerm(velocity=𝐮,scheme=WENO5())
terms  = (term1,)
# msh = Mesh.GenericMesh{2,Float64}()
# LevelSetMethods.meshgen!(msh,ls)
# plot(msh)

# solve
integrator = RK2()
eq = LevelSetEquation(;terms,integrator,levelset=deepcopy(ϕ),t=0,cfl=0.5,boundary_condition=bc)
anim = @animate for n ∈ 0:100
    tf = 0.05*n
    integrate!(eq,tf)
    plot(eq)
end
gif(anim, "test.gif")
