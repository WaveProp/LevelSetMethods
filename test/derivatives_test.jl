using Test
using LinearAlgebra
using StaticArrays

import WavePropBase as WPB
import LevelSetMethods as LSM

@testset "Uniform mesh" begin
    rec   = WPB.HyperRectangle((-2.,-2.),(2.,2.))
    ϕ    = LSM.CartesianGridFunction(rec;meshsize=0.1) do (x,y)
        x^2 + y^2 - 1
    end
    ∇ϕ   = LSM.CartesianGridFunction(rec;meshsize=0.1) do (x,y)
        SVector(2x,2y)
    end
    h = step(ϕ)
    # first derivative
    I  = CartesianIndex(9,7)
    for op in (LSM.D⁺,LSM.D⁻,LSM.D⁰,LSM.weno5⁻,LSM.weno5⁺)
        for dir in 1:2
            ee =  abs(∇ϕ[I][dir] - op(ϕ,I,dir))
            @test ee < 5*h[dir]
        end
    end
    # second derivative, same direction
    for dir in 1:2
        @test abs( 2 - LSM.D2⁰(ϕ,I,dir)) < 5*h[dir]
        @test abs( 2 - LSM.D2(ϕ,I,(dir,dir))) < 5*h[dir]
    end
    # second derivative, different directions
    for op in (LSM.D2,)
        for dims in ((1,2),(2,1))
            @test abs(op(ϕ,I,dims)) < 5*h[dims[1]]*h[dims[2]]
        end
    end
end
