using SafeTestsets
using LevelSetMethods

@safetestset "Bernstein polynomial" begin include("bernsteinpolynomial_test.jl") end
@safetestset "Grid functions" begin include("gridfunction_test.jl") end
@safetestset "Derivatives" begin include("derivatives_test.jl") end
# TODO: add tests for boundary conditions
@safetestset "Level set" begin include("levelset_test.jl") end
@safetestset "Mesh generation" begin include("meshgen_test.jl") end
