"""
    struct MeshField{V,M}
    
A field described by its discrete values on a mesh.
"""
struct MeshField{V,M}
    vals::V
    mesh::M
end

"""
    LevelSet

Alias for [`MeshField`](@ref).
"""
const LevelSet{V,M} = MeshField{V,M}

function MeshField(f::Function,m)
    vals = map(f,m)
    MeshField(vals,m)
end    

# getters
mesh(ϕ::MeshField) = ϕ.mesh
Base.values(ϕ::MeshField) = ϕ.vals

# geometric dimension
dimension(f::MeshField) = dimension(mesh(f))

# overload base methods for convenience
Base.getindex(ϕ::MeshField,I...) = getindex(values(ϕ),I...)
Base.setindex!(ϕ::MeshField,vals,I...) = setindex!(values(ϕ),vals,I...)
Base.size(ϕ::MeshField) = size(values(ϕ))
Base.eltype(ϕ::MeshField) = eltype(values(ϕ))
Base.zero(ϕ::MeshField) = MeshField(zero(values(ϕ)),mesh(ϕ))

# apply b.c. on a MeshField by rewritting the "border nodes"
function applybc!(mf::MeshField,bctype)
    # periodic, only one layer of ghost nodes   
    ϕ = values(mf) 
    if bctype == :periodic1 
        # @. ϕ[1,:]   = ϕ[end-1,:]
        @views copy!(ϕ[1,:],ϕ[end-1,:])
        @views copy!(ϕ[end,:],ϕ[2,:])
        @views copy!(ϕ[:,1],ϕ[:,end-1])
        @views copy!(ϕ[:,end],ϕ[:,2])
    elseif bctype == :neumann1
        @views ϕ[1,:]   = ϕ[2,:]
        ϕ[end,:] = ϕ[end-1,:]
        ϕ[:,1]   = ϕ[:,2]
        ϕ[:,end] = ϕ[:,end-1]
    else 
        error("uknown boundary condition $bctype")    
    end    
    return ϕ
end    

# some helper functions
function _increment_index(I::CartesianIndex,dim::Integer)
    N = length(I)    
    @assert 1 ≤ dim ≤ length(I)
    return I + CartesianIndex(ntuple(i-> i==dim,N))
end    

function _decrement_index(I::CartesianIndex,dim::Integer)
    N = length(I)    
    @assert 1 ≤ dim ≤ length(I)
    return I + CartesianIndex(ntuple(i-> -(i==dim),N))
end    

# recipes for Plots
@recipe function f(ϕ::MeshField)
    N = dimension(ϕ)    
    if N == 2 # 2d contour plot
        seriestype := :contour
        levels --> [0,]
        m = mesh(ϕ)    
        return xgrid(m),ygrid(m),values(ϕ)
    else
        notimplemented()
    end        
end    
