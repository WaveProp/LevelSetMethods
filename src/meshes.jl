"""
    abstract type AbstractMesh{N,T}
    
An abstract mesh structure in dimension `N` with primite data of type `T`. 
"""
abstract type AbstractMesh{N,T} end

struct CartesianGrid{N,T} <: AbstractMesh{N,T}
    grid1d::NTuple{N,Vector{T}}
end    

# allow for arguments to be passed as e.g. (x,y) instead of tuple and promote
# type if needed
CartesianGrid(args...) = CartesianGrid(promote(args...))

Base.size(g::CartesianGrid) = map(length,g.grid1d)

meshsize(g::CartesianGrid)  = map(x->x[2]-x[1],g.grid1d)


