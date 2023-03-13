"""
    abstract type BoundaryCondition

Singleton types used in [`applybc!`](@ref) for dispatch purposes.

Boundary conditions are imposed by a *ghost node* method, meaning that meshes
are composed of *internal nodes*, and *ghost nodes* which sole purpose is to
reinforce a given boundary condition.
"""
abstract type BoundaryCondition end

"""
    num_ghost_points(bc::BoundaryCondition)

The number of ghost-points on which to apply `bc`.
"""
function num_ghost_points end

"""
    applybc!(ϕ,bc::BoundaryCondition,[dim])

Overwrite the ghost nodes of `ϕ` do impose the boundary condition `bc` on the
dimension `dim`. If no `dim` argument is passed, the boundary condition is
applied to all dimensions of `ϕ`.
"""
function applybc! end

function applybc!(ϕ::CartesianGridFunction,bc::BoundaryCondition)
    N = ambient_dimension(ϕ)
    for dim in 1:N
        applybc!(ϕ,bc,dim)
    end
    return ϕ
end

"""
    interior_indices(ϕ,bc::BoundaryCondition)

Indices to iterate over the interior nodes of `ϕ`.
"""
function interior_indices end

function interior_indices(ϕ,bc::BoundaryCondition)
    P  = num_ghost_points(bc)
    sz = size(ϕ)
    N  = length(sz)
    # range of nodes in each dimenion. Excludes the first P and last P nodes
    I  = ntuple(N) do dim
        P+1:sz[dim]-P
    end
    return CartesianIndices(I)
end

"""
    struct PeriodicBC{P}

A periodic boundary condition with `P` periodic nodes.
"""
struct PeriodicBC{N} <: BoundaryCondition end

PeriodicBC(n::Int) = PeriodicBC{n}()

num_ghost_points(::PeriodicBC{N}) where {N} = N

function applybc!(ϕ::CartesianGridFunction,bc::PeriodicBC{P},dir) where {P}
    A = vals(ϕ)
    I⁻r,I⁺r = _index_read(A,bc,dir)
    I⁻w,I⁺w = _index_write(A,bc,dir)
    @views copy!(A[I⁻w],A[I⁺r])
    @views copy!(A[I⁺w],A[I⁻r])
    return ϕ
end


"""
    struct NeumannBC{P}

A high-order homogeneous Neumann boundary condition with `P` ghost nodes.
"""
struct NeumannBC{N} <: BoundaryCondition end

NeumannBC(n::Int) = NeumannBC{n}()

num_ghost_points(::NeumannBC{N}) where {N} = N

function applybc!(ϕ::CartesianGridFunction,bc::NeumannBC{P},dir) where {P}
    A = vals(ϕ)
    I⁻r,I⁺r = _index_read(A,bc,dir)
    I⁻w,I⁺w = _index_write(A,bc,dir)
    @views copy!(A[I⁻w],A[I⁻r])
    @views copy!(A[I⁺w],A[I⁺r])
    return ϕ
end

## Helper functions

function _index_write(ϕ,bc::BoundaryCondition,dir)
    P = num_ghost_points(bc)
    sz = size(ϕ)
    N  = length(sz)
    Il = ntuple(N) do dim
        if dim == dir
            1:P
        else
            1:sz[dim]
        end
    end
    Ir = ntuple(N) do dim
        if dim == dir
            (sz[dim]-P+1):sz[dim]
        else
            1:sz[dim]
        end
    end
    return CartesianIndices(Il), CartesianIndices(Ir)
end

function _index_read(ϕ,bc::BoundaryCondition,dir)
    P  = num_ghost_points(bc)
    sz = size(ϕ)
    N  = length(sz)
    Il = ntuple(N) do dim
        if dim == dir
            (P+1):2P
        else
            1:sz[dim]
        end
    end
    Ir = ntuple(N) do dim
        if dim == dir
            (sz[dim]-2P+1):(sz[dim]-P)
        else
            1:sz[dim]
        end
    end
    return CartesianIndices(Il), CartesianIndices(Ir)
end
