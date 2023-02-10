using Test
using LevelSetMethods
import WavePropBase as WPB

import LevelSetMethods as LSM

WPB.clear_entities!()
box = WPB.HyperRectangle((-1,-1),(1,1))
# create a circle of radius one
f = x -> x[1]^2 + x[2]^2 - 1
Ω = LSM.LevelSet(f,box,-1)
@test WPB.ambient_dimension(Ω) == 2
@test WPB.geometric_dimension(Ω) == 2

@test WPB.ambient_dimension(WPB.boundary(Ω)) == 2
@test WPB.geometric_dimension(WPB.boundary(Ω)) == 1

fh  = LSM.CartesianGridFunction(f,box; meshsize=0.1,order=4)
Ωh   = LSM.LevelSet(fh)
@test WPB.ambient_dimension(Ωh) == 2
@test WPB.geometric_dimension(Ωh) == 2

@test WPB.ambient_dimension(WPB.boundary(Ωh)) == 2
@test WPB.geometric_dimension(WPB.boundary(Ωh)) == 1
