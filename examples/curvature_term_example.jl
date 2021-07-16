using Test
using LevelSetMethods
using WavePropBase
using LinearAlgebra
using Plots

hx,hy = 0.05,0.05
domain = HyperRectangle((-1,-1),(1,1))
M      = UniformCartesianMesh(domain;step=(hx,hy))
bc   = PeriodicBC(1)
ϕ    = LevelSet(M,bc) do (x,y)
    Inf
end
add_circle!(ϕ,SVector(0.5,0.0),0.25)
add_circle!(ϕ,SVector(-0.5,0.0),0.25)
add_rectangle!(ϕ,SVector(0.0,0.0),SVector(1.0,0.1))
b     = NodeField(M) do (x,y)
    -min(hx,hy)
end
term  = CurvatureTerm(b)
integrator = ForwardEuler()
eq = LevelSetEquation(;terms=(term,),integrator,state=ϕ,t=0)

dt = 0.01
anim = @animate for n ∈ 0:100
    tf = dt*n
    integrate!(eq,tf)
    fig = plot(eq,linecolor=:black,linestyle = :solid)
    # fig = heatmap(eq,colorbar=true)
end
gif(anim, "test.gif")
