using Test
using StaticArrays
import WavePropBase as WPB
import LevelSetMethods as LSM

U = WPB.HyperRectangle((0.,-1.,-1.), (2.,2.,3.))
L = WPB.low_corner(U); H = WPB.high_corner(U)
a = zeros(3,3,3)
a[1,1,1] = -1.; a[3,1,1] = 1.; a[1,3,1] = 2.; a[1,1,3] = 3.
ϕ = LSM.power2bernstein(a, U) # x1^2 + 2x2^2 + 3x3^2 - 1
N = 5 # number of tests of each category

@testset "Algebra tests" begin
    @testset "Evaluation tests" begin
        f = (x) -> x[1]^2 + 2x[2]^2 + 3x[3]^2 - 1
        for _ in 1:N
            x = @SVector rand(3)
            x = x .* (H - L) .+ L
            @test @inferred(ϕ(x)) ≈ f(x)
        end
    end

    @testset "Derivation tests" begin
        ∇ϕ = LSM.partials(ϕ)
        ∇f = SVector(x->2x[1], x->4x[2], x->6x[3])
        for _ in 1:N
            x = @SVector rand(3)
            x = x .* (H - L) .+ L
            for i in 1:3
                @test @inferred(∇ϕ[i](x)) ≈ ∇f[i](x)
            end
        end
    end

    @testset "Splitting tests" begin
        c = WPB.center(U)
        for d in 1:3
            ϕ1, ϕ2 = split(ϕ, d)
            for _ in 1:N
                x = WPB.svector(i->i == d ? c[i] : rand()*(H[i]-L[i])+L[i], 3)
                x̂ = deleteat(x, d)
                @test LSM.upper_restrict(ϕ1, d)(x̂) ≈ LSM.lower_restrict(ϕ2, d)(x̂) ≈ ϕ(x)
            end
        end
    end
end
