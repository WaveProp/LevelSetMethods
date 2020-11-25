using Test
using LevelSetMethods
using LinearAlgebra
using Plots

nx,ny = 100,100
x     = LinRange(-1,1,nx)
y     = LinRange(-1,1,ny)
hx,hy = step(x),step(y)
grid = CartesianGrid(x,y)
ϕ    = LevelSet(grid) do (x,y)
    1.0
end
add_circle!(ϕ,SVector(0.5,0.0),0.25)
add_circle!(ϕ,SVector(-0.5,0.0),0.25)
add_rectangle!(ϕ,SVector(0.0,0.0),SVector(1.0,0.1))
plot(ϕ)
v     = MeshField(grid) do (x,y)
    -1.0
end
b     = MeshField(x->-0.5,grid)
term1  = NormalAdvectionTerm(v)
terms = (term1,)
# need at least two ghost cells for second order scheme for normal advection
bc    = PeriodicBC(2)
pgap = 1
# integrator = ForwardEuler(cfl=0.5,buffer=similar(ϕ))
integrator = RK2(cfl=0.5,buffers=(similar(ϕ),similar(ϕ)))
t = 0
anim = @animate for n ∈ 0:20
    ϕ, t = LevelSetMethods.evolve!(ϕ,integrator,terms,bc,t)
    if n%pgap == 0
        plot(ϕ,title="t=$t")
    end
end

gif(anim, "test.gif", fps = 15)
