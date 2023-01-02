function marchingcubes(ls::CartesianLevelSet)
    @assert ambient_dimension(ls) == 3 "marching cubes only works in 3d"
    f   = levelset_function(ls)
    msh = f.vals_mesh
    x,y,z = collect.(grids(msh))
    mc = MC(f.vals,Int;x,y,z)
    march(mc)
    T = Triangle3D{Float64}
    msh = GenericMesh{3,Float64}()
    append!(msh.nodes,mc.vertices)
    connectivity = collect(reinterpret(reshape,Int64,mc.triangles))
    push!(msh.elements,T=>connectivity)
    push!(msh.ent2tags,ls=>Dict(T=>[i for i in 1:size(connectivity,2)]))
    return msh
end
