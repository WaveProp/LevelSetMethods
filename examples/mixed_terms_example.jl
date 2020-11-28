using Test
using LevelSetMethods
using LinearAlgebra
using Plots

nx,ny = 100,100
x     = LinRange(-1,1,nx)
y     = LinRange(-1,1,ny)
hx,hy = step(x),step(y)
grid = CartesianGrid(x,y)
bc    = PeriodicBC(3)
ϕ    = LevelSet(grid,bc) do (x,y)
    1.0
end
add_circle!(ϕ,SVector(0.5,0.0),0.25)
add_circle!(ϕ,SVector(-0.5,0.0),0.25)
add_rectangle!(ϕ,SVector(0.0,0.0),SVector(1.0,0.1))
plot(ϕ)
v     = MeshField(grid) do (x,y)
    -0.1
end
𝐮     = MeshField(grid) do (x,y)
    SVector(-y,x)
end   
b     = MeshField(grid) do (x,y)
    -min(hx,hy)
end   
term1  = NormalMotionTerm(v)
term2  = AdvectionTerm(velocity=𝐮,scheme=WENO5())
# term2  = AdvectionTerm(velocity=𝐮,scheme=Upwind())
term3  = CurvatureTerm(b)
terms = (term1,term2,term3)
buf     = zero(ϕ)
integrator = ForwardEuler(0.4)
eq = LevelSetEquation(;terms,integrator,state=ϕ,t=0,buffer=buf)

dt = 0.01
anim = @animate for n ∈ 0:80
    tf = dt*n    
    integrate!(eq,tf)    
    plot(eq,linecolor=:black)    
end
gif(anim, "test.gif")
