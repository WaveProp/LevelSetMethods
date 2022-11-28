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
    TensorLagInterp,
    NodeIterator,
    LagrangeSquare,
    LagrangeTriangle,
    ParametricElement,
    ReferenceHyperCube,
    SType,
    PolynomialSpace,
    monomial_basis,
    mesh,
    ent2tags,
    domain,
    grids,
    nodes,
    boundary,
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

# static functionality
include("gridfunction.jl")
include("levelset.jl")
include("linearization.jl")
include("meshgen.jl")

# dynamic functionality
# include("boundaryconditions.jl")
# include("derivatives.jl")
# include("timestepping.jl")
# include("levelsetequation.jl")

export
    HyperRectangle,
    LevelSetEquation,
    LevelSet,
    CartesianLevelSet,
    CartesianGridFunction,
    PeriodicBC,
    AdvectionTerm,
    CurvatureTerm,
    Upwind,
    WENO5,
    ForwardEuler,
    RK2,
    integrate!,
    meshgen,
    boundary

end # module
