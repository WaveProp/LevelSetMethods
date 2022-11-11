using Test
using LevelSetMethods
using LinearAlgebra
using StaticArrays
import WavePropBase as WPB

using LevelSetMethods: D⁺, D⁻, D⁰, D2⁰,D2, weno5⁻, weno5⁺

@testset "Uniform mesh" begin
    nx,ny = 100,50
    rec   = WPB.HyperRectangle((-2.,-2.),(2.,2.))
    grid  = WPB.UniformCartesianMesh(rec,(nx,ny))
    h    = step(grid)
    ϕ    = CartesianGridFunction(grid) do (x,y)
        x^2 + y^2 - 1
    end
    ∇ϕ   = CartesianGridFunction(grid) do (x,y)
        SVector(2x,2y)
    end
    # first derivative
    I  = CartesianIndex(9,7)
    for op in (D⁺,D⁻,D⁰,weno5⁻,weno5⁺)
        for dir in 1:2
            ee =  abs(∇ϕ[I][dir] - op(ϕ,I,dir))
            @test ee < 5*h[dir]
        end
    end
    # second derivative, same direction
    for dir in 1:2
        @test abs( 2 - D2⁰(ϕ,I,dir)) < 5*h[dir]
        @test abs( 2 - D2(ϕ,I,(dir,dir))) < 5*h[dir]
    end
    # second derivative, different directions
    for op in (D2,)
        for dims in ((1,2),(2,1))
            @test abs(op(ϕ,I,dims)) < 5*h[dims[1]]*h[dims[2]]
        end
    end
end
