"""
    struct NodeField{V,M}

A field described by its discrete `vals::V` on a `mesh::M`. A boundary condition
`bc` may aditionally be specified, in which case the nodes of the `NodeField`
can be [`interior_nodes`](@ref) or `boundary_nodes`.
"""
struct NodeField{V,M,B}
    vals::V
    mesh::M
    bc::B
end

# getters
mesh(ϕ::NodeField) = ϕ.mesh
Base.values(ϕ::NodeField) = ϕ.vals
grid1d(mf::NodeField,args...) = grid1d(mesh(mf),args...)

has_boundary_condition(mf::NodeField) = (mf.bc !== nothing)
boundary_condition(mf) = mf.bc

function NodeField(f::Function,m)
    iter = NodeIterator(m)
    vals = map(f,iter)
    NodeField(vals,m,nothing)
end

# geometric dimension
ambient_dimension(f::NodeField) = ambient_dimension(mesh(f))

function Base.step(f::NodeField{<:Any,<:UniformCartesianMesh})
    msh = mesh(f)
    step(msh)
end
Base.step(f::NodeField,i::Int) = step(f)[i]

# overload base methods for convenience
Base.getindex(ϕ::NodeField,I...) = getindex(values(ϕ),I...)
Base.setindex!(ϕ::NodeField,vals,I...) = setindex!(values(ϕ),vals,I...)
Base.size(ϕ::NodeField) = size(values(ϕ))
Base.eltype(ϕ::NodeField) = eltype(values(ϕ))
Base.zero(ϕ::NodeField) = NodeField(zero(values(ϕ)),mesh(ϕ),boundary_condition(ϕ))
Base.similar(ϕ::NodeField) = NodeField(similar(values(ϕ)),mesh(ϕ),boundary_condition(ϕ))

"""
    LevelSet

Alias for [`NodeField`](@ref) with a boundary condition `bc::BoundaryCondition`.
"""
const LevelSet{V,M,B<:BoundaryCondition} = NodeField{V,M,B}

function LevelSet(f::Function,m,bc::BoundaryCondition=PeriodicBC(0))
    iter = NodeIterator(m)
    vals = map(f,iter)
    ϕ    = NodeField(vals,m,bc)
    # applybc!(ϕ)
    return ϕ
end

"""
    applybc!(ϕ::LevelSet)

Overwrite the boundary nodes of `ϕ` to impose the boundary conditio `ϕ.bc`.
"""
applybc!(ϕ::LevelSet) = applybc!(ϕ,boundary_condition(ϕ))

function interior_indices(ϕ::LevelSet)
    iter = mesh(ϕ) |> NodeIterator
    interior_indices(iter,boundary_condition(ϕ))
end

# helps to obtain classical shapes's signed distance function
function circle_signed_distance(m,center,r)
    rsq = r*r
    return map(x->sum((x.-center).^2)-rsq,m)
end

function rectangle_signed_distance(m,center,size)
    sized2 = 0.5*size
    return map(x->maximum(abs.(x.-center)-sized2),m)
end

# helpers to add geometric shapes on the a level set
function add_circle!(ϕ::NodeField,center,r)
    iter = mesh(ϕ) |> NodeIterator
    circle = circle_signed_distance(iter,center,r)
    union!(values(ϕ),circle)
end
function remove_circle!(ϕ::NodeField, center, r)
    iter = mesh(ϕ) |> NodeIterator
    circle = circle_signed_distance(iter,center,r)
    difference!(values(ϕ),circle)
end

function add_rectangle!(ϕ::NodeField, center, size)
    iter = mesh(ϕ) |> NodeIterator
    rectangle = rectangle_signed_distance(iter,center,size)
    union!(values(ϕ),rectangle)
end
function remove_rectangle!(ϕ::NodeField, center, size)
    iter = mesh(ϕ) |> NodeIterator
    rectangle = rectangle_signed_distance(iter,center,size)
    difference!(values(ϕ),rectangle)
end

# helpers to merge or make the difference between two level set functions
@inline function Base.union!(ϕ1,ϕ2)
    @. ϕ1 = min(ϕ1,ϕ2)
end
@inline function difference!(ϕ1,ϕ2)
    @. ϕ1 = -min(-ϕ1,ϕ2)
end

# recipes for Plots
@recipe function f(ϕ::NodeField)
    N = ambient_dimension(ϕ)
    if N == 2 # 2d contour plot
        seriestype --> :contour
        levels --> [0,]
        aspect_ratio --> :equal
        colorbar --> false
        # seriescolor --> :black
        m = mesh(ϕ)
        # Note: the values of ϕ need be transposed because contour expects the
        # matrix to have rows representing the x values and columns expecting
        # the y value.
        xgrid = grids(m)[1]
        ygrid = grids(m)[2]
        return xgrid,ygrid,transpose(values(ϕ))
    else
        notimplemented()
    end
end
