module LevelSetMethods

using LinearAlgebra
using StaticArrays
using RecipesBase
using Base.Threads: @threads, @spawn

using WavePropBase
using WavePropBase.Geometry
using WavePropBase.Mesh
using WavePropBase.Utils
WavePropBase.@import_interface

export
    UniformCartesianMesh,
    HyperRectangle,
    SVector,
    NodeField,
    LevelSet,
    PeriodicBC,
    AdvectionTerm,
    CurvatureTerm,
    NormalMotionTerm,
    ReinitializationTerm,
    compute_terms,
    add_circle!,
    remove_circle!,
    add_rectangle!,
    remove_rectangle!,
    ForwardEuler,
    RK2,
    RKLM2,
    Upwind,
    WENO5,
    LevelSetEquation,
    NodeIterator,
    ElementIterator,
    integrate!

include("boundaryconditions.jl")
include("nodefield.jl")
include("derivatives.jl")
include("levelsetterms.jl")
include("timestepping.jl")
include("levelsetequation.jl")

end # module
