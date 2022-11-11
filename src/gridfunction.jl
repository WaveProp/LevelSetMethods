"""
    struct CartesianGridFunction{N,T,V}

A function defined by discrete values on a the nodes of a `UniformCartesianMesh`.
"""
struct CartesianGridFunction{N,T,V} <: AbstractArray{V,N}
    vals::Array{V,N}
    mesh::UniformCartesianMesh{N,T}
end

mesh(f::CartesianGridFunction) = f.mesh
vals(f::CartesianGridFunction) = f.vals
Base.step(f::CartesianGridFunction, args...) = step(mesh(f), args...)
ambient_dimension(f::CartesianGridFunction{N}) where {N} = N

# AbstractArray interface
Base.size(f::CartesianGridFunction) = size(vals(f))
Base.getindex(f::CartesianGridFunction, args...) = getindex(vals(f), args...)
Base.setindex!(f::CartesianGridFunction, args...) = setindex!(vals(f), args...)
Base.eltype(f::CartesianGridFunction) = eltype(vals(f))

domain(f::CartesianGridFunction) = domain(mesh(f))

function CartesianGridFunction(f::Function, msh::UniformCartesianMesh)
    vals = [f(x) for x in NodeIterator(msh)]
    CartesianGridFunction(vals, msh)
end

function (f::CartesianGridFunction{N})(x::SVector{N}) where {N}
    @assert N == 2
    m = mesh(f)
    T = eltype(f)
    I = element_index_for_point(x, m)
    p = interpolant(f, I, LagrangeSquare{4,T})
    # p     = interpolant(f,I,LagrangeTriangle{6,T})
    return p(x)
end

function interpolant(f::CartesianGridFunction, I, ::Type{<:LagrangeSquare{4}})
    # TODO: properly implement other interpolants
    m = mesh(f)
    els = ElementIterator(m)
    rec = els[I]
    I1, I2, I3, I4 = I, I + CartesianIndex(1, 0), I + CartesianIndex(1, 1), I + CartesianIndex(0, 1)
    vals = SVector(f[I1], f[I2], f[I3], f[I4])
    p̂ = LagrangeSquare(vals)
    p = (x) -> begin
        x̂ = (x - low_corner(rec)) ./ width(rec) # normalized coordinates
        p̂(x̂)
    end
    return p
end

function interpolant(f::CartesianGridFunction, I, ::Type{<:LagrangeTriangle{3}})
    m = mesh(f)
    els = ElementIterator(m)
    rec = els[I]
    I1, I2, I3 = I, I + CartesianIndex(1, 0), I + CartesianIndex(0, 1)
    vals = SVector(f[I1], f[I2], f[I3])
    p̂ = LagrangeTriangle(vals)
    p = (x) -> begin
        x̂ = (x - low_corner(rec)) ./ width(rec) # normalized coordinates in [0,1] × [0,1]
        p̂(x̂)
    end
    return p
end

function interpolant(f::CartesianGridFunction, I, ::Type{<:LagrangeTriangle{6}})
    m = mesh(f)
    els = ElementIterator(m)
    rec = els[I]
    I1, I2, I3 = I, I + CartesianIndex(2, 0), I + CartesianIndex(0, 2)
    I4, I5, I6 = I + CartesianIndex(1, 0), I + CartesianIndex(1, 1), I + CartesianIndex(0, 1)
    vals = SVector(f[I1], f[I2], f[I3], f[I4], f[I5], f[I6])
    p̂ = LagrangeTriangle(vals)
    p = (x) -> begin
        x̂ = (x - low_corner(rec)) ./ width(rec) ./ 2 # normalized coordinates in [0,1] × [0,1]
        p̂(x̂)
    end
    return p
end
