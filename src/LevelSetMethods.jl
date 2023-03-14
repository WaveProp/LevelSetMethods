module LevelSetMethods

import WavePropBase as WPB

using LinearAlgebra
using Roots
using Contour
using MarchingCubes
using IntervalArithmetic
using ForwardDiff
using StaticArrays
using RecipesBase
using Printf

import WavePropBase:
    AbstractEntity,
    HyperRectangle,
    UniformCartesianMesh,
    NystromMesh,
    GenericMesh,
    ElementIterator,
    TensorLagInterp,
    GaussLegendre,
    QuadratureNode,
    NodeIterator,
    LagrangeSquare,
    LagrangeTriangle,
    ParametricElement,
    ReferenceHyperCube,
    SType,
    PolynomialSpace,
    monomial_basis,
    mesh,
    degree,
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
    meshgen,
    meshgen!,
    integrate

include("utils.jl")
include("bernsteinpolynomials.jl")
include("gridfunction.jl")
include("derivatives.jl")
include("boundaryconditions.jl")
include("levelset.jl")
include("meshgen.jl")
include("quadgen.jl")
include("timestepping.jl")
include("levelsetequation.jl")
include("plotsIO.jl")

export
    LevelSet,
    CartesianLevelSet,
    meshgen


# export
#     NystromMesh,
#     LevelSetEquation,
#     LevelSet,
#     CartesianLevelSet,
#     CartesianGridFunction,
#     PeriodicBC,
#     AdvectionTerm,
#     CurvatureTerm,
#     Upwind,
#     WENO5,
#     ForwardEuler,
#     RK2,
#     integrate,
#     meshgen,
#     boundary,
#     power2bernstein

end # module
