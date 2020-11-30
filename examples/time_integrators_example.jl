using Test
using LevelSetMethods
using LinearAlgebra
using Plots

nx,ny = 50,50
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
term2  = AdvectionTerm(velocity=𝐮)
term3  = CurvatureTerm(b)
terms = (term1,term2,term3)
b = (zero(ϕ),zero(ϕ))
integrator = ForwardEuler()
eq = LevelSetEquation(;terms,integrator,state=ϕ,t=0,buffer=b[1])
integrator = RK2()
eq2 = LevelSetEquation(;terms,integrator,state=deepcopy(ϕ),t=0,buffer=deepcopy(b))
integrator = LevelSetMethods.RKLM2()
eq3 = LevelSetEquation(;terms,integrator,state=deepcopy(ϕ),t=0,buffer=deepcopy(b[1]))

dt = 0.01
anim = @animate for n ∈ 0:50
    tf = dt*n    
    integrate!(eq,tf)    
    integrate!(eq2,tf)    
    integrate!(eq3,tf)    
    fig = plot(eq,linecolor=:black,linestyle = :dash)
    plot!(fig,eq2,linecolor=:blue,linestyle = :dot)
    plot!(fig,eq3,linecolor=:red,linestyle = :dashdot)
    # pplot(eq,eq2)
end
gif(anim, "test.gif")


