"""
    abstract type LevelSetTerm
    
A typical term in a level-set evolution equation. 

These are functor-like structures, callable as `term(ϕ::LevelSet)`, and
returning the desired approximation with the same type as `ϕ`. Additionally, a
buffer may be passed in the functor call as `term(buffer,ϕ) == buffer
+ term(ϕ)`.
"""
abstract type LevelSetTerm end

"""
    struct AdvectionTerm{V,M} <: LevelSetTerm

Level-set advection term such that `(adv::AdvectionTerm)(ϕ) ≈  𝐯 ⋅ ∇ϕ`. The velocity
field `𝐯` is represented as a [`MeshField`](@ref).
"""
@Base.kwdef struct AdvectionTerm{V,M} <: LevelSetTerm
    velocity::MeshField{V,M}
    bc::Symbol = :periodic1
    scheme::Symbol = :upwind
end

function (adv::AdvectionTerm)(buffer::MeshField,ϕ::MeshField;scheme=:upwind)
    @assert mesh(adv.velocity) == mesh(ϕ) == mesh(buffer)
    𝐮 = adv.velocity # velocity field
    bc = adv.bc
    if adv.scheme == :upwind
        advect_upwind!(buffer,ϕ,𝐮,bc) # an advection scheme. The signature is fixed.
    else
        notimplemented()    
    end
    return buffer
end    
(adv::AdvectionTerm)(ϕ::MeshField;scheme=advect_upwind!) = adv(zero(ϕ),ϕ;scheme)

function advect_upwind!(buffer::MeshField,ϕ::MeshField,𝐮::MeshField,bc)
    @assert mesh(ϕ) === mesh(𝐮)
    grid = mesh(ϕ)
    N    = dimension(ϕ)
    sz   = size(grid)
    h    = meshsize(grid)
    applybc!(ϕ,bc)
    for I in CartesianIndices(grid)
        # check for border cases
        isborder = any(1:N) do dim
            i = I[dim]
            i == 1 || i == sz[dim]
        end
        isborder && continue
        # for each dimension, compute the upwind derivative and multiply by
        # velocity and add to buffer
        for dim in 1:N
            v = 𝐮[I][dim]
            if v > 0
                Im = _decrement_index(I,dim)    
                buffer[I] += v*(ϕ[I] - ϕ[Im]) / h[dim]    
            else
                Ip = _increment_index(I,dim)    
                buffer[I] += v*(ϕ[Ip] - ϕ[I]) / h[dim]    
            end    
        end
        # buffer[I] += sum(1:N) do dim   
        #     v = 𝐮[I][dim]
        #     if v > 0
        #         Im = _decrement_index(I,dim)    
        #         v*(ϕ[I] - ϕ[Im]) / h[dim]    
        #     else
        #         Ip = _increment_index(I,dim)    
        #         v*(ϕ[Ip] - ϕ[I]) / h[dim]    
        #     end    
        # end
    end        
    return buffer
end

"""
    struct CurvatureTerm{V,M} <: LevelSetTerm

Level-set curvature term such that `(curv::CurvatureTerm)(ϕ) ≈  b κ |∇ϕ|`, where
`κ = ∇ ⋅ (∇ϕ / |∇ϕ|)`. The scalar field `b` is represented as a [`MeshField`](@ref).
"""
struct CurvatureTerm{V,M} <: LevelSetTerm
    b::MeshField{V,M}
end    






