using SafeTestsets
using LevelSetMethods

@safetestset "Grid functions" begin include("gridfunction_test.jl") end

@safetestset "Level set" begin include("levelset_test.jl") end

# @safetestset "Level set" begin include("linearization_test.jl") end

# @safetestset "Derivatives" begin include("derivatives_test.jl") end

@safetestset "Mesh generation" begin include("meshgen_test.jl") end
