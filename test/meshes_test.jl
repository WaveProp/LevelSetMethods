using Test
using LevelSetMethods

@testset "Basic ops" begin
    nx,ny = 100,50
    x     = LinRange(-1,1,nx)
    y     = LinRange(0,3,ny)
    grid  = UniformCartesianMesh(x,y)
    iter  = NodeIterator(grid)
    @test size(iter) === (length(x), length(y))
    @test step(grid)[1] ≈ step(x)
    @test step(grid)[2] ≈ step(y)
    @test all(iter[6,2] .== (x[6],y[2]))
    @test length(CartesianIndices(iter)) == length(x)*length(y)
end
