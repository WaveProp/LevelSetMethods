#=
    Implementation of methods required to generate conforming elements from a
    level-set by recursive subdivision and one-dimensional root-finding.
=#

function meshgen(ls::LevelSet; maxdepth=20, maxslope=10, meshsize, order=1)
    meshgen(CartesianLevelSet(ls; meshsize, order=order);maxdepth, maxslope)
end

function meshgen(ls::CartesianLevelSet; maxdepth=20, maxslope=10)
    N = ambient_dimension(ls)
    msh = GenericMesh{N,Float64}()
    p   = (;maxdepth,maxslope)
    return meshgen!(msh, ls, p)
end

function meshgen!(msh::GenericMesh, ls::CartesianLevelSet, p::NamedTuple)
    N = ambient_dimension(ls)
    V = SVector{N,Float64}
    D = ReferenceHyperCube{Int(geometric_dimension(ls))} # reference domain
    s = levelset_sign(ls)
    edict = msh.elements
    e2t = Dict{DataType,Vector{Int}}()
    # loop over the C∞ interpolants of `ls`
    for f in bernstein_interpolants(ls)
        root = MultiBernsteinCell([f], [s])
        surf = s == 0
        level = 0
        maps  = _meshgen(root, surf, level, p)
        # sort elements by type for the given entity
        for τ in maps
            el = ParametricElement{D,V}(τ)
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
        union!(roots, find_zeros(ψ, L, U))
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

function _meshgen(Ω::MultiBernsteinCell{N,T}, surf, level,p) where {N,T}
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
                y = find_zero(ψₖ, (lk, rk))
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
