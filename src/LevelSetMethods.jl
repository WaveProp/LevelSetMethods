module LevelSetMethods

using LinearAlgebra
using StaticArrays
using RecipesBase

export 
    CartesianGrid, 
    meshsize,
    SVector,
    MeshField,
    LevelSet,
    PeriodicBC,
    AdvectionTerm,
    CurvatureTerm,
    compute_terms
    
include("meshes.jl")
include("boundaryconditions.jl")
include("meshfield.jl")
include("derivatives.jl")
include("levelsetterms.jl")
include("timestepping.jl")

end # module
