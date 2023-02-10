#=
    Implementation of methods required to generate conforming elements from a
    level-set by recursive subdivision and one-dimensional root-finding.
=#

struct LevelSetElement{L,F,D,T} <: WPB.AbstractElement{D,T}
    ϕ::L # level-set function
    parametrization::F # explicit parametrization mapping D into the element
    function LevelSetElement{D,T}(ϕ::L, f::F) where {L,F,D,T}
        return new{L,F,D,T}(ϕ, f)
    end
end

WPB.ambient_dimension(el::LevelSetElement{L,F,D,T}) where {L,F,D,T} = length(T)
function WPB.geometric_dimension(el::LevelSetElement{L,F,D,T}) where {L,F,D,T}
    return geometric_dimension(D)
end

(el::LevelSetElement)(u::SVector) = el.parametrization(u)
(el::LevelSetElement)(u::Tuple) = el.parametrization(SVector(u))

function WPB.normal(el::LevelSetElement, u)
    x = el(u)
    ∇ϕ = ForwardDiff.gradient(el.ϕ, x)
    return normalize(∇ϕ)
end

function WPB.curvature(el::LevelSetElement, u)
    x = el(u)
    ∇ϕ = ForwardDiff.gradient(el.ϕ, x)
    Hp = ForwardDiff.hessian(el.ϕ, x)
    Np = norm(∇ϕ)
    return tr(Hp) / Np - dot(∇ϕ, Hp, ∇ϕ) / Np^3
end

function WPB.meshgen(ls::LevelSet; kwargs...)
    N = ambient_dimension(ls)
    T = Float64
    msh = GenericMesh{N,T}()
    return meshgen!(msh, ls; kwargs...)
end

function WPB.meshgen!(msh, ls::LevelSet; maxdepth=20, maxslope=10)
    ϕ = levelset_function(ls)
    @assert ϕ isa CartesianGridFunction "the `levelset_function` must be a `CartesianGridFunction`"
    s = levelset_sign(ls)
    p = (; maxdepth, maxslope, s)
    return _meshgen!(msh, ls, ϕ, p)
end

function _meshgen!(msh::GenericMesh{N,T}, ls::LevelSet, ϕ::CartesianGridFunction{N,T},
                   p::NamedTuple) where {N,T}
    V = SVector{N,T}
    D = ReferenceHyperCube{Int(geometric_dimension(ls))} # reference domain
    edict = msh.elements
    e2t = Dict{DataType,Vector{Int}}()
    # loop over the C∞ interpolants of `ls`
    for f in bernstein_interpolants(ϕ)
        root = MultiBernsteinCell([f], [p.s])
        surf = p.s == 0
        level = 0
        maps = _meshgen(root, surf, level, p)
        # sort elements by type for the given entity
        for τ in maps
            el = LevelSetElement{D,V}(f, τ)
            E = typeof(el)
            els = get!(edict, E, Vector{E}())
            tags = get!(e2t, E, Int[])
            push!(els, el)
            push!(tags, length(els))
        end
    end
    # store maps in the mesh
    haskey(ent2tags(msh), ls) && (@warn "remeshed $ls")
    ent2tags(msh)[ls] = e2t
    return msh
end

function dim1mesh(Ψ::Vector, signs::Vector{<:Integer}, L, U)
    roots = [L, U]
    for ψ in Ψ
        y = my_find_zero(ψ, (L, U))
        isnothing(y) || union!(roots, y)
        # union!(roots, find_zeros(ψ, L, U))
    end
    sort!(roots)
    Maps = Vector{Function}()
    for (l, r) in zip(roots[1:(end - 1)], roots[2:end])
        val = [ψ((l + r) / 2.0) for ψ in Ψ]
        if all((val .* signs) .≥ 0)
            τ(x) = x * (r - l) .+ l
            push!(Maps, τ)
        end
    end
    return Maps
end

function HDmesh(Ψ::Vector, signs::Vector{<:Integer}, L, U, τ, k, D)
    roots = [_ -> L, _ -> U]
    x₀ = svector(_ -> 0.5, D - 1)
    x̂₀ = τ(x₀)
    for ψ in Ψ
        if ψ(insert(x̂₀, k, L)) * ψ(insert(x̂₀, k, U)) < 0
            push!(roots, x -> find_zero(y -> ψ(insert(τ(x), k, y)), (L, U)))
        end
    end
    sort!(roots; lt=(l, r) -> l(x₀) < r(x₀))
    Maps = Vector{Function}()
    for (l, r) in zip(roots[1:(end - 1)], roots[2:end])
        val = [ψ(insert(x̂₀, k, (l(x₀) + r(x₀)) / 2.0)) for ψ in Ψ]
        if all((val .* signs) .≥ 0)
            function t(x)
                x₁ = deleteat(x, k)
                x̂₁ = τ(x₁)
                return insert(x̂₁, k, x[k] * (r(x₁) - l(x₁)) + l(x₁))
            end
            push!(Maps, t)
        end
    end
    return Maps
end

function _meshgen(Ω::MultiBernsteinCell{N,T}, surf, level, p) where {N,T}
    rec = Ω.rec
    xc = center(rec)
    Ψ = Ω.Ψ
    ∇Ψ = Ω.∇Ψ
    signs = Ω.signs
    D = ambient_dimension(rec)
    @assert !surf || (D > 1)
    ctype = cell_type(Ω)
    # check for limiting cases of empty or whole cells
    if ctype == empty_cell
        return Vector{Function}()
    elseif ctype == whole_cell
        if surf
            return Vector{Function}()
        else
            return [x -> x .* (high_corner(rec) .- low_corner(rec)) .+ low_corner(rec)]
        end
    end
    # return if max level reached
    if level ≥ p.maxdepth
        if !isempty(Ψ)
            @show first(Ψ).domain
        end
        @warn "Maximum depth $(p.maxdepth) reached: resorting to low-order approximation"
        if surf
            return [_ -> xc]
        else
            return [_ -> xc]
        end
    end

    # base case
    if D == 1
        return dim1mesh(Ψ, signs, low_corner(rec)[1], high_corner(rec)[1])
    end

    # find a heigh direction such that all of ∇Ψ are (provably) bounded away
    # from zero.
    bnds = map(∇Ψ) do ∇ψ
        return ntuple(d -> bound(∇ψ[d], rec), D)
    end
    isvalid = ntuple(D) do dim
        all(bnds) do bnd
            return (prod(bnd[dim]) > 0) &&
                   (sum(bd -> maximum(abs, bd), bnd) / minimum(abs, bnd[dim]) <
                    p.maxslope)
        end
    end
    if !any(isvalid) # no valid direction so split
        Ω1, Ω2 = split(Ω)
        Maps1 = _meshgen(Ω1, surf, level + 1, p)
        Maps2 = _meshgen(Ω2, surf, level + 1, p)
        test = Vector{Function}()
        return append!(test, Maps1, Maps2)
    end

    # If there is a valid direction, we go down on it. Choose the direction which
    # is the least steep overall by maximizing the minimum of the derivative on
    # direction k over all functions
    ∇Ψc = map(∇Ψ) do ∇ψ
        ntuple(D) do d
            ∇ψc = abs.(∇ψ[d](xc))
            return ∇ψc / norm(∇ψc)
        end
    end
    k = argmax(1:D) do dim
        if isvalid[dim]
            minimum(∇ψc -> abs(∇ψc[dim]), ∇Ψc)
        else
            -Inf
        end
    end
    Ω̃ = restrict(Ω, k, surf) # the D-1 dimensional domain
    maps = _meshgen(Ω̃, false, level, p)
    Maps = Vector{Function}()
    for t in maps
        if surf
            # if we get here there should be only one Ψ
            @assert length(Ψ) == 1
            ψ = first(Ψ)
            lk, rk = low_corner(rec)[k], high_corner(rec)[k]
            function τ(x)
                x̂ = t(x)
                ψₖ(y) = ψ(insert(x̂, k, y))
                y = my_find_zero(ψₖ, (lk, rk))
                isnothing(y) && error()
                return insert(x̂, k, y)
            end
            push!(Maps, τ)
        else
            lk, rk = low_corner(rec)[k], high_corner(rec)[k]
            M = HDmesh(Ψ, signs, lk, rk, t, k, D)
            append!(Maps, M)
        end
    end
    ###########################################
    return Maps
end

## Marching squares
function marchingsquares(ls::LevelSet)
    msh = GenericMesh{2,Float64}()
    return marchingsquares!(msh, ls)
end

function marchingsquares!(msh::GenericMesh{N,T}, ls::LevelSet) where {N,T}
    @assert ambient_dimension(ls) == 2 "marching squares requires two-dimensional level set"
    @assert levelset_sign(ls) == 0 "marching squares can only mesh zero level sets: got s = $(levelset_sign(ϕ))"
    ϕ = levelset_function(ls)
    @assert ϕ isa CartesianGridFunction "marching squares requires the level-set function to be a `CartesianGridFunction`"
    return _marchingsquares!(msh, ls, ϕ)
end

function _marchingsquares!(msh, ls::LevelSet, ϕ::CartesianGridFunction)
    x, y = WPB.grids(vals_mesh(ϕ))
    z = ϕ.vals
    # elements are lines in 2d with two points
    E = WPB.LagrangeLine{2,WPB.Point2D}
    enttags = Int[]
    msh.ent2tags[ls] = Dict(E => enttags)
    el2tag = get!(msh.elements, E, Matrix{Int}(undef, 2, 0))
    elidx = size(el2tag, 2)
    el2tag = vec(el2tag)
    cl = contour(x, y, z, 0)
    for line in lines(cl)
        pts = line.vertices
        closed = pts[1] == pts[end]
        npts = length(pts)
        istart = length(msh.nodes) + 1
        for n in 1:(npts - 1)
            push!(msh.nodes, pts[n])
            i = length(msh.nodes)
            push!(el2tag, i, i + 1)
            elidx += 1
            push!(enttags, elidx)
        end
        if closed
            # make last node index point to first
            el2tag[end] = istart
        else
            push!(msh.nodes, pts[end])
        end
    end
    el2tag = reshape(el2tag, 2, :)
    msh.elements[E] = el2tag
    return msh
end

## Marching cubes
function marchingcubes(ls::LevelSet)
    @assert ambient_dimension(ls) == 3 && geometric_dimension(ls) == 2 "marching cubes only works on three-dimensional surfaces"
    N, T = 3, Float64
    msh = GenericMesh{N,T}()
    return marchingcubes!(msh, ls)
end

function marchingcubes!(msh::GenericMesh, ls)
    ϕ = levelset_function(ls)
    msg = "marching cubes requires the level-set function to be a `CartesianGridFunction`"
    @assert ϕ isa CartesianGridFunction msg
    x, y, z = collect.(grids(vals_mesh(ϕ)))
    mc = MC(vals(ϕ), Int; x, y, z)
    march(mc)
    T = WPB.Triangle3D{Float64}
    nstart = length(msh.nodes)
    append!(msh.nodes, mc.vertices)
    connectivity = nstart .+ collect(reinterpret(reshape, Int64, mc.triangles))
    if haskey(msh.elements, T)
        msh.elements[T] = hcat(msh.elements[T], connectivity)
    else
        msh.elements[T] = connectivity
    end
    push!(msh.ent2tags, ls => Dict(T => [i for i in 1:size(connectivity, 2)]))
    return msh
end
