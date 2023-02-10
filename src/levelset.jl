"""
    struct LevelSet <: AbstractEntity

A geometrical entity implicitly represented through a function-like object `f :
ℝᵈ → ℝ`. The underlying region is given by `Ω = {𝐱 ∈ U : s*f(𝐱)>0}`, for `s ∈
{+1,-1}` (i.e. a volume), and by `Γ = {𝐱 ∈ U : f(𝐱)=0}` (i.e. a surface) if
`s=0`, where `U` is a `HyperRectangle` used to provide a `bounding_box` for the
domain.

```jldoctest
using LevelSetMethods
box = HyperRectangle((-2,-2),(2,2))
# create a disk of radius one
f = x -> x[1]^2 + x[2]^2 - 1
Ω = LevelSet(f,box)
```
"""
struct LevelSet <: AbstractEntity
    dim::UInt8 # geometrical dimension
    tag::Int
    f # maps ℝᵈ → ℝ, where d is the ambient dimension
    s::Int
    U::HyperRectangle
    boundary::Union{LevelSet,Nothing}
    function LevelSet(d::Integer, t::Integer, f, s, U, bnd)
        @assert (s == 0 || s == -1 || s == 1) "s must be either 0, -1 or 1"
        ent = new(d, t, f, s, U, bnd)
        # every entity gets added to a global variable ENTITIES so that we can
        # ensure the (d,t) pair is a UUID for an entity, and to easily retrieve
        # different entities.
        global_add_entity!(ent)
        return ent
    end
end

bounding_box(ent::LevelSet) = ent.U
ambient_dimension(ent::LevelSet) = ambient_dimension(bounding_box(ent))
levelset_function(ent::LevelSet) = ent.f
levelset_sign(ent::LevelSet) = ent.s

function LevelSet(f, U::HyperRectangle, s=-1)
    if s == 0
        dim = ambient_dimension(U) - 1
        tag = new_tag(dim)
        bnd = nothing
        return LevelSet(dim, tag, f, s, U, bnd)
    else
        dim = ambient_dimension(U)
        tag = new_tag(dim)
        bnd = LevelSet(f, U, 0)
        return LevelSet(dim, tag, f, s, U, bnd)
    end
end
LevelSet(f::CartesianGridFunction,s::Int=-1) = LevelSet(f, domain(f), s)

function Base.union(ϕ1::LevelSet, ϕ2::LevelSet)
    @assert ϕ1.s == ϕ2.s "the two level set entities must have the same sign"
    box = bounding_box(ϕ1) ∪ bounding_box(ϕ2)
    f = (x) -> min(ϕ1.f(x), ϕ2.f(x))
    return LevelSet(f, ϕ1.s, box)
end

# TODO: implement something similar to "simple_shapes" of parametric surfaces
