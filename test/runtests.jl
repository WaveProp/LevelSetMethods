using SafeTestsets
using LevelSetMethods

@safetestset "Grid functions" begin include("gridfunction_test.jl") end

@safetestset "Level set" begin include("levelsetentity_test.jl") end

@safetestset "Derivatives" begin include("derivatives_test.jl") end
