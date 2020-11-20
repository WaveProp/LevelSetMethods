module LevelSetMethods

using LinearAlgebra
using StaticArrays

export 
    CartesianGrid, 
    meshsize,
    upwind,
    evolve!
    
include("meshes.jl")
include("levelsets.jl")
# include("derivatives.jl")
# include("timestepping.jl")

end # module
