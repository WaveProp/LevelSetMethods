using Test
using LinearAlgebra
import WavePropBase as WPB
import LevelSetMethods as LSM

@testset "Marching squares" begin
    # circle
    box = WPB.HyperRectangle((-2.0,-2),(2,2))
    f = x -> x[1]^2 + x[2]^2 - 1
    ϕ = LSM.CartesianGridFunction(f,box;meshsize=0.1)
    Ω = LSM.LevelSet(ϕ,box)
    Γ = WPB.boundary(Ω)
    M = LSM.marchingsquares(Γ)
    Q = WPB.NystromMesh(M;qorder=1)
    # since the approximation is piecewise linear, can't ask for too much from the perimeter
    @test abs(WPB.integrate(x->1,Q) - 2π) < 1e-2
    # two (open) curves
    f = x -> abs(x[2]) - 0.5
    ϕ = LSM.CartesianGridFunction(f,box;meshsize=0.4)
    Ω = LSM.LevelSet(ϕ,box)
    Γ = WPB.boundary(Ω)
    M = LSM.marchingsquares(Γ)
    Q = WPB.NystromMesh(M;qorder=1)
    # the domain are two straigh-lines, so the linear approximation is exact
    @test WPB.integrate(x->1,Q) ≈ 8
end

@testset "Marching cubes" begin
    # sphere
    box = WPB.HyperRectangle((-2.,-2,-2.),(2,2,2))
    f = x -> norm(x) - 1
    ϕ = LSM.CartesianGridFunction(f,box;meshsize=0.11)
    Ω1 = LSM.LevelSet(ϕ,box)
    Γ1 = WPB.boundary(Ω1)
    M = LSM.marchingcubes(Γ1)
    #WPB.vtk_mesh_file(M,"sphere") |> vtk_save
    Q = WPB.NystromMesh(M;qorder=1)
    # since the approximation is piecewise linear, can't ask for too much from the perimeter
    @test abs(WPB.integrate(x->1,Q) - 4π)/4π < 1e-2
    # add another sphere to the mesh
    box = WPB.HyperRectangle((0,-2,-2.),(4,2,2))
    f = x -> norm(x .- (2,0,0)) - 1
    ϕ = LSM.CartesianGridFunction(f,box;meshsize=0.11)
    Ω2 = LSM.LevelSet(ϕ,box)
    Γ2 = WPB.boundary(Ω2)
    LSM.marchingcubes!(M,Γ2)
    Q = WPB.NystromMesh(M;qorder=1)
    @test abs(WPB.integrate(x->1,Q) - 8π)/8π < 1e-2
    # two (open) curves
    f = x -> abs(x[2]) - 0.5
    ϕ = LSM.CartesianGridFunction(f,box;meshsize=0.4)
    Ω = LSM.LevelSet(ϕ,box,0)
    M = LSM.marchingcubes(Ω)
    Q = WPB.NystromMesh(M;qorder=1)
    # WPB.vtk_mesh_file(M,"planes") |> vtk_save
    # the domain are two straigh-lines, so the linear approximation is exact
    @test WPB.integrate(x->1,Q) ≈ 32
end

@testset "meshgen" begin
    WPB.clear_entities!()
    box = WPB.HyperRectangle((-2.0,-2),(2,2))
    # create a circle of radius one
    f = x -> x[1]^2 + x[2]^2 - 1
    ϕ = LSM.CartesianGridFunction(f,box;meshsize=0.1,order=2)
    Ω = LSM.LevelSet(ϕ,box,-1)
    M = LSM.meshgen(Ω)
    Q = WPB.NystromMesh(M;qorder=10)

    @test WPB.ambient_dimension(Ω) == 2
    @test WPB.geometric_dimension(Ω) == 2
    @test WPB.integrate(x->1,Q) ≈ π

    Γ = WPB.boundary(Ω)
    M = LSM.meshgen(Γ)
    Q = WPB.NystromMesh(M;qorder=20)

    @test WPB.ambient_dimension(WPB.boundary(Ω)) == 2
    @test WPB.geometric_dimension(WPB.boundary(Ω)) == 1
    @test_broken WPB.integrate(x->1,Q) ≈ 2π

    # 3d
    # FIXME: the 3d tests are too slow so the tolerance has been dropped to 1
    # percent for testing areas/volumes. Once we work on the performance, fix these
    # tests to a more stringent tolerance like in 2d.
    # WPB.clear_entities!()
    # l = 1.0
    # box = WPB.HyperRectangle((-l,-l,-l),(l,l,l))
    # # create a circle of radius one
    # ϕ = LSM.CartesianGridFunction(box;meshsize=0.1,order=2) do x
    #     sqrt(x[1]^2 + x[2]^2 + x[3]^2) - 1
    # end
    # Ω = LevelSet(ϕ)
    # M = meshgen(Ω)
    # Q = WPB.NystromMesh(M;qorder=3)

    # @test WPB.ambient_dimension(Ω) == 3
    # @test WPB.geometric_dimension(Ω) == 3
    # @test abs(WPB.integrate(x->1,Q) - 4/3*π) < 1e-2

    # Γ = WPB.boundary(Ω)
    # M = meshgen(Γ;meshsize=0.5,order=3)
    # Q = WPB.NystromMesh(M;qorder=4)

    # @test WPB.ambient_dimension(WPB.boundary(Ω)) == 3
    # @test WPB.geometric_dimension(WPB.boundary(Ω)) == 2
    # @test abs(WPB.integrate(x->1,Q) - 4π) < 1e-2

end
