using Test
using LevelSetMethods

@testset "Basic ops" begin
    hx,hy = 0.1, 0.2
    x = collect(-1:hx:1)
    y = collect(0:hy:3)
    grid = CartesianGrid(x,y)
    @test size(grid) === (length(x), length(y))
    @test meshsize(grid)[1] ≈ hx
    @test meshsize(grid)[2] ≈ hy
    @test grid[6,2] == (x[6],y[2])
    @test length(CartesianIndices(grid)) == length(x)*length(y)
end

