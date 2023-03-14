using Test
using LinearAlgebra
using StaticArrays
using Plots

import WavePropBase as WPB
import LevelSetMethods as LSM

WPB.clear_entities!()
nx,ny = 200,200
rec   = WPB.HyperRectangle((-3,-3),(3.,3.))
ϕ     = LSM.CartesianGridFunction((x) -> sqrt((x[1]-1.5)^2 + x[2]^2) - 1, rec; meshsize=0.1)
bc    = LSM.PeriodicBC(3)
Γ     = LSM.LevelSet(ϕ,0)
nodes = WPB.NodeIterator(LSM.vals_mesh(ϕ))
u⃗     = [2π*SVector(-x[2],x[1]) for x in nodes]
term1 = LSM.AdvectionTerm(velocity=u⃗,scheme=LSM.WENO5())
terms = (term1,)
integrator = LSM.RK2()
eq = LSM.LevelSetEquation(;terms,integrator,levelset=Γ,t=0,cfl=0.5,boundary_condition=bc)

# solve
anim = @animate for n ∈ 0:100
    tf = 0.05*n
    LSM.integrate!(eq,tf)
    plot(eq)
end
gif(anim, "test.gif")
