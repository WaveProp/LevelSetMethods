"""
    LevelSet <: AbstractEntity

A geometrical entity implicitly represented through a function `f : ℝᵈ → ℝ`. The
underlying region is given by `Ω = {𝐱 ∈ U : s*f(𝐱)>0}`, for `s ∈ {+1,-1}` (i.e. a
volume), and by `Γ = {𝐱 ∈ U : f(𝐱)=0}` (i.e. a surface) if `s=0`, where `U` is
a `HyperRectangle`.
"""
struct LevelSet <: AbstractEntity
    dim::UInt8 # geometrical dimension
    tag::Int
    f::Function
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

function LevelSet(f, s::Int, U::HyperRectangle)
    if s == 0
        dim = ambient_dimension(U) - 1
        tag = new_tag(dim)
        bnd = nothing
        return LevelSet(dim, tag, f, s, U, bnd)
    else
        dim = ambient_dimension(U)
        tag = new_tag(dim)
        bnd = LevelSet(f, 0, U)
        return LevelSet(dim, tag, f, s, U, bnd)
    end
end

function Base.union(ϕ1::LevelSet, ϕ2::LevelSet)
    @assert ϕ1.s == ϕ2.s "the two level set entities must have the same sign"
    box = bounding_box(ϕ1) ∪ bounding_box(ϕ2)
    f = (x) -> min(ϕ1.f(x), ϕ2.f(x))
    return LevelSet(f, ϕ1.s, box)
end

# TODO: implement other set operatos

struct DiscreteLevelSet <: AbstractEntity
    dim::UInt8 # geometrical dimension
    tag::Int
    f::CartesianGridFunction
    s::Int
    boundary::Union{DiscreteLevelSet,Nothing}
    function DiscreteLevelSet(d::Integer, t::Integer, f, s, bnd)
        @assert (s == 0 || s == -1 || s == 1) "s must be either 0, -1 or 1"
        ent = new(d, t, f, s, bnd)
        # every entity gets added to a global variable ENTITIES so that we can
        # ensure the (d,t) pair is a UUID for an entity, and to easily retrieve
        # different entities.
        global_add_entity!(ent)
        return ent
    end
end

function DiscreteLevelSet(f::CartesianGridFunction, s::Int)
    U = domain(f)
    if s == 0
        dim = ambient_dimension(U) - 1
        tag = new_tag(dim)
        bnd = nothing
        return DiscreteLevelSet(dim, tag, f, s, bnd)
    else
        dim = ambient_dimension(U)
        tag = new_tag(dim)
        bnd = DiscreteLevelSet(f, 0)
        return DiscreteLevelSet(dim, tag, f, s, bnd)
    end
end

function DiscreteLevelSet(f::Function,msh,s=-1)
    g = CartesianGridFunction(f,msh)
    return DiscreteLevelSet(g,s)
end

(ϕ::DiscreteLevelSet)(x) = ϕ.f(x)

vals(ϕ::DiscreteLevelSet) = vals(ϕ.f)
mesh(ϕ::DiscreteLevelSet) = mesh(ϕ.f)
bounding_box(ent::DiscreteLevelSet) = domain(mesh(ent))
ambient_dimension(ent::DiscreteLevelSet) = ambient_dimension(bounding_box(ent))
Base.step(ϕ::DiscreteLevelSet,args...) = step(ϕ.f,args...)
Base.size(ϕ::DiscreteLevelSet) = size(ϕ.f)
Base.getindex(ϕ::DiscreteLevelSet, args...) = getindex(ϕ.f, args...)
Base.setindex!(ϕ::DiscreteLevelSet, args...) = setindex!(ϕ.f, args...)

# TODO: implement something similar to "simple_shapes" of parametric surfaces


# recipes for Plots
@recipe function f(ls::DiscreteLevelSet)
    ϕ = ls.f
    s = ls.s
    N = ambient_dimension(ϕ)
    if N == 2 # 2d contour plot
        if s == 0
            seriestype --> :contour
        else
            seriestype --> :contourf
        end
        levels --> [0,0]
        aspect_ratio --> :equal
        colorbar --> false
        # seriescolor --> :black
        m = mesh(ϕ)
        # Note: the vals of ϕ need be transposed because contour expects the
        # matrix to have rows representing the x vals and columns expecting
        # the y value.
        xgrid = grids(m,1)
        ygrid = grids(m,2)
        return xgrid,ygrid,-transpose(vals(ϕ))
    else
        notimplemented()
    end
end
