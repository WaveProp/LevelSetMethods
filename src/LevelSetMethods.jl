module LevelSetMethods

using LinearAlgebra
using StaticArrays
using RecipesBase
using Base.Threads: @threads, @spawn

import WavePropBase:
    AbstractEntity,
    HyperRectangle,
    UniformCartesianMesh,
    ElementIterator,
    NodeIterator,
    LagrangeSquare,
    LagrangeTriangle,
    mesh,
    domain,
    grids,
    element_index_for_point,
    vals,
    ambient_dimension,
    new_tag,
    low_corner,
    width,
    global_add_entity!,
    increment_index,
    decrement_index

include("boundaryconditions.jl")
include("gridfunction.jl")
include("levelset.jl")
include("derivatives.jl")
include("timestepping.jl")
include("levelsetequation.jl")
# include("meshgen.jl")

export
    LevelSetEquation,
    LevelSet,
    DiscreteLevelSet,
    CartesianGridFunction,
    PeriodicBC,
    AdvectionTerm,
    Upwind,
    WENO5,
    ForwardEuler,
    RK2,
    integrate!

end # module
