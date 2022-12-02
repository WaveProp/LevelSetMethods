using Test
using LevelSetMethods
import WavePropBase as WPB

WPB.clear_entities!()
box = WPB.HyperRectangle((-2.0,-2),(2,2))
# create a circle of radius one
f = x -> x[1]^2 + x[2]^2 - 1
Ω  = LevelSet(f,box,-1)
M = meshgen(Ω; meshsize=0.1, order=2)
Q = WPB.NystromMesh(M;qorder=10)

@test WPB.ambient_dimension(Ω) == 2
@test WPB.geometric_dimension(Ω) == 2
@test WPB.integrate(x->1,Q) ≈ π

Ωc = CartesianLevelSet(Ω;meshsize=0.1,order=3)
M  = meshgen(Ωc)
Q  = WPB.NystromMesh(M;qorder=20)

@test WPB.ambient_dimension(Ω) == 2
@test WPB.geometric_dimension(Ω) == 2
@test WPB.integrate(x->1,Q) ≈ π

Γ = WPB.boundary(Ω)
M = meshgen(Γ;meshsize=0.1,order=2)
Q = WPB.NystromMesh(M;qorder=20)

@test WPB.ambient_dimension(WPB.boundary(Ω)) == 2
@test WPB.geometric_dimension(WPB.boundary(Ω)) == 1
@test WPB.integrate(x->1,Q) ≈ 2π

Γ = WPB.boundary(Ωc)
M = meshgen(Γ)
Q = WPB.NystromMesh(M;qorder=20)

@test WPB.ambient_dimension(WPB.boundary(Ω)) == 2
@test WPB.geometric_dimension(WPB.boundary(Ω)) == 1
@test WPB.integrate(x->1,Q) ≈ 2π

# 3d
# FIXME: the 3d tests are too slow so the tolerance has been dropped to 1
# percent for testing areas/volumes. Once we work on the performance, fix these
# tests to a more stringent tolerance like in 2d.
WPB.clear_entities!()
l = 1.0
box = WPB.HyperRectangle((-l,-l,-l),(l,l,l))
# create a circle of radius one
f = x -> sqrt(x[1]^2 + x[2]^2 + x[3]^2) - 1
Ω = LevelSet(f,box,-1)
M = meshgen(Ω;meshsize=0.5,order=3)
Q = WPB.NystromMesh(M;qorder=3)

@test WPB.ambient_dimension(Ω) == 3
@test WPB.geometric_dimension(Ω) == 3
@test abs(WPB.integrate(x->1,Q) - 4/3*π) < 1e-2

Γ = WPB.boundary(Ω)
M = meshgen(Γ;meshsize=0.5,order=3)
Q = WPB.NystromMesh(M;qorder=4)

@test WPB.ambient_dimension(WPB.boundary(Ω)) == 3
@test WPB.geometric_dimension(WPB.boundary(Ω)) == 2
@test abs(WPB.integrate(x->1,Q) - 4π) < 1e-2
