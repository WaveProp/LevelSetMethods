using Test
using LevelSetMethods
import WavePropBase as WPB

WPB.clear_entities!()
box = WPB.HyperRectangle((-1,-1),(1,1))
# create a circle of radius one
f = x -> x[1]^2 + x[2]^2 - 1
Ω = LevelSet(f,-1,box)
@test WPB.ambient_dimension(Ω) == 2
@test WPB.geometric_dimension(Ω) == 2

@test WPB.ambient_dimension(WPB.boundary(Ω)) == 2
@test WPB.geometric_dimension(WPB.boundary(Ω)) == 1
