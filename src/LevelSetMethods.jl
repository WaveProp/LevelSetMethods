module LevelSetMethods

using LinearAlgebra
using StaticArrays
using RecipesBase

export 
    CartesianGrid, 
    meshsize,
    upwind,
    evolve!,
    SVector,
    MeshField,
    LevelSet,
    AdvectionTerm
    
include("meshes.jl")
include("meshfield.jl")
include("levelsetterms.jl")

end # module
