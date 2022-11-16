module LevelSetMethods

using LinearAlgebra
using StaticArrays
using Roots
using RecipesBase
using Base.Threads: @threads, @spawn

import WavePropBase:
    AbstractEntity,
    HyperRectangle,
    UniformCartesianMesh,
    GenericMesh,
    ElementIterator,
    NodeIterator,
    LagrangeSquare,
    LagrangeTriangle,
    ParametricElement,
    ReferenceHyperCube,
    mesh,
    ent2tags,
    domain,
    grids,
    element_index_for_point,
    vals,
    ambient_dimension,
    geometric_dimension,
    gradient,
    new_tag,
    svector,
    center,
    half_width,
    low_corner,
    high_corner,
    section,
    width,
    global_add_entity!,
    increment_index,
    decrement_index,
    meshgen

include("boundaryconditions.jl")
include("gridfunction.jl")
include("levelset.jl")
include("linearization.jl")
include("meshgen.jl")
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
    integrate!,
    meshgen

end # module
