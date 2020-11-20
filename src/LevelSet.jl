module LevelSet

using LinearAlgebra
using StaticArrays

export 
    CartesianGrid, 
    upwind,
    evolve!
    
include("cartesiangrid.jl")
include("derivatives.jl")
include("timestepping.jl")

end # module
