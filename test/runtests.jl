using SafeTestsets
using LevelSetMethods

@safetestset "Derivatives tests" begin include("derivatives_test.jl") end

@safetestset "Meshes tests" begin include("meshes_test.jl") end
