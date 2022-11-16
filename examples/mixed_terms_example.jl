using Test
using LevelSetMethods
using LinearAlgebra
using Plots

hx,hy = 0.01,0.02
domain = HyperRectangle((-1,-1),(1,1))
M      = UniformCartesianMesh(domain;step=(hx,hy))
bc   = PeriodicBC(3)
ls    = LevelSet(M) do (x,y)
    Inf
end
add_circle!(ls.ϕ,SVector(0.5,0.0),0.25)
add_circle!(ls.ϕ,SVector(-0.5,0.0),0.25)
add_rectangle!(ls.ϕ,SVector(0.0,0.0),SVector(1.0,0.1))
v     = NodeField(M) do (x,y)
    -0.1
end
u     = NodeField(M) do (x,y)
    SVector(-y,x)
end
b     = NodeField(M) do (x,y)
    -min(hx,hy)
end
term1  = NormalMotionTerm(v)
term2  = AdvectionTerm(velocity=u)
term3  = CurvatureTerm(b)
terms = (term1,term2,term3)
terms = (term1,)
integrator = ForwardEuler()
eq = LevelSetEquation(;terms,integrator,state=ϕ,t=0,boundary_condition=bc)

dt = 0.01
anim = @animate for n ∈ 0:100
    tf = dt*n
    integrate!(eq,tf)
    fig = plot(eq,linecolor=:black,linestyle = :solid)
end
gif(anim, "test.gif")
