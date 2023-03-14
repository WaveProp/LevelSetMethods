using Test
using LinearAlgebra
using StaticArrays
using Plots

import WavePropBase as WPB
import LevelSetMethods as LSM

WPB.clear_entities!()
nx,ny = 200,200
rec   = WPB.HyperRectangle((-3,-3),(3.,3.))
ϕ     = LSM.CartesianGridFunction((x) -> sqrt(x[1]^2 + x[2]^2) - 1, rec; meshsize=0.1)
bc    = LSM.NeumannBC(3)
Γ     = LSM.LevelSet(ϕ,0)
term1 = LSM.ReinitializationTerm()
terms = (term1,)
integrator = LSM.ForwardEuler()
eq = LSM.LevelSetEquation(;terms,integrator,levelset=Γ,t=0,cfl=0.05,boundary_condition=bc)

fig = heatmap(eq,colorbar=true)
plot!(eq,color=:black,colorbar=true)

LSM.integrate!(eq,1)

@show LSM.df_deviation(eq)

# solve
anim = @animate for n ∈ 0:100
    tf = 0.05*n
    @show LSM.df_deviation(eq)
    LSM.integrate!(eq,tf)
    plot(eq)
end
gif(anim, "test.gif")
