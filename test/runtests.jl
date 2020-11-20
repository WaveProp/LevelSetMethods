using SafeTestsets
using LevelSetMethods

@safetestset "Meshes" begin include("meshes_test.jl") end

@safetestset "Advection term" begin include("advect_test.jl") end
