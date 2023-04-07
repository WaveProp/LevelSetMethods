using Test
using LevelSetMethods
using LevelSetMethods: CartesianGridFunction, interpolant, merge_interpolant
import WavePropBase as WPB

for dim in 2:3
    box = WPB.HyperRectangle(ntuple(i->-i,dim),ntuple(i->2,dim))
    meshsize = 0.1
    # a bilinear function, interpolation should be exact
    f1 =  x -> 3 + 2*x[1] - x[2] + 0.5*x[1]*x[2]
    f2 =  x -> 3 + 2*x[1] - x[2] + 0.5*x[1]*x[2] + 0.5*x[1]^2 + 0.5*x[dim]^2 + prod(xd->xd^2,x)
    f3 =  x -> 3 + 2*x[1] - x[2] + 0.5*x[1]*x[2] + 0.5*x[1]^2 + 0.5*x[2]^2 - 0.1*x[dim]^3 + prod(xd->xd^3,x)
    for (k,f) in enumerate((f1,f2,f3))
        for order in 1:3
            F = CartesianGridFunction(f,box;meshsize,order)
            I = CartesianIndex(ntuple(i->i+1,dim))
            # p = interpolant(F, I)
            # U = WPB.domain(p)
            # xt = U(ntuple(i->rand(),dim))
            # @testset "test function of order $k, dim $dim, interp order $order" begin
            #     order < k ? (@test p(xt) ≉ f(xt)) : (@test p(xt) ≈ f(xt))
            # end
            for J in LevelSetMethods.neighbors(I, size(mesh(F)))
                p = merge_interpolant(F, I, J)
                U = WPB.domain(p)
                xt = U(ntuple(i->rand(),dim))
                @testset "test merge of order $k, dim $dim, interp order $order" begin
                    order < k ? (@test p(xt) ≉ f(xt)) : (@test p(xt) ≈ f(xt))
                end
            end
           
        end
    end
end
