using Test
using LevelSetMethods
using LevelSetMethods: CartesianGridFunction
import WavePropBase as WPB

WPB.clear_entities!()
box = WPB.HyperRectangle((-2,-2),(2,2))

msh = WPB.UniformCartesianMesh(box,(10,10))

els   = WPB.ElementIterator(msh)
nodes = WPB.NodeIterator(msh)

# a bilinear function, interpolation should be exact
f = x -> 3 + 2*x[1] - x[2] + 0.5*x[1]*x[2]

F = CartesianGridFunction(f,msh)
for _ in 1:5
    x₀ = rand(WPB.Point2D)
    @test F(x₀) ≈ f(x₀)
end

f = x -> x[1]^2 + x[2]^2 - 1

F = CartesianGridFunction(f,msh)
