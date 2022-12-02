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

@enum CellType empty_cell whole_cell cut_cell

struct MultiBernsteinCell{N,T}
    Ψ::Vector{BernsteinPolynomial{N,T}}
    ∇Ψ::Vector{SVector{N,BernsteinPolynomial{N,T}}}
    signs::Vector{Int}
    rec::HyperRectangle{N,T}
    celltype::CellType
    function MultiBernsteinCell(Ψ::Vector{BernsteinPolynomial{N,T}}, signs) where {N,T}
        @assert length(signs) == length(Ψ)
        rec = domain(first(Ψ))
        @assert all(ψ -> domain(ψ) == rec, Ψ)
        ∇Ψ = map(partials, Ψ)
        ctype = _prune!(Ψ, ∇Ψ, rec, signs)
        return new{N,T}(Ψ, ∇Ψ, signs, rec::HyperRectangle{N,T}, ctype)
    end
end

function MultiBernsteinCell(ψ::BernsteinPolynomial, s::Integer; kwargs...)
    return MultiBernsteinCell([ψ], [s]; kwargs...)
end

cell_type(Ω::MultiBernsteinCell) = Ω.celltype

"""
    prune!(Ω)

Prune the functions specifying the domain `Ω` and return the `CellType` of the domain.
"""
function prune!(Ω::MultiBernsteinCell)
    return _prune!(Ω.Ψ, Ω.∇Ψ, Ω.rec, Ω.signs)
end

function _prune!(Ψ, ∇Ψ, rec, signs)
    delInd = Vector{Int}()
    for (i, ψ) in enumerate(Ψ)
        si = signs[i]
        t = cell_type(ψ, si, rec)
        if t == whole_cell
            # intersection is the whole rec, so ψ can be prune
            # @info "Whole cell"
            append!(delInd, i)
        elseif t == empty_cell
            # intersection is empty, return immediately
            return empty_cell
        end
    end
    deleteat!(signs, delInd)
    deleteat!(Ψ, delInd)
    deleteat!(∇Ψ, delInd)
    isempty(Ψ) && return whole_cell
    return cut_cell
end

function cell_type(ψ::BernsteinPolynomial, s, rec)
    l, u = bound(ψ)
    ψc = ψ(center(rec))
    l * u ≥ 0 || (return cut_cell)
    if s * ψc ≥ 0
        # intersection is the whole rec
        return whole_cell
    else
        # intersection is empty, return immediately
        return empty_cell
    end
end

function lower_restrict(ψ, rec, k)
    a = low_corner(rec)[k]
    return x -> ψ(insert(x, k, a))
end

function upper_restrict(ψ, rec, k)
    a = high_corner(rec)[k]
    return x -> ψ(insert(x, k, a))
end

function lower_restrict_grad(∇ψ::SVector{N}, rec, k) where {N}
    a = low_corner(rec)[k]
    ∇ψ′ = deleteat(∇ψ, k)
    (x) -> ∇ψ′(insert(x, k, a))
    return svector(d -> (x) -> ∇ψ′[d](insert(x, k, a)), N - 1)
end

function upper_restrict_grad(∇ψ::SVector{N}, rec, k) where {N}
    a = high_corner(rec)[k]
    ∇ψ′ = deleteat(∇ψ, k)
    return svector(d -> (x) -> ∇ψ′[d](insert(x, k, a)), N - 1)
end

function Base.split(Ω::MultiBernsteinCell)
    Ψ = Ω.Ψ
    signs = Ω.signs
    rec = Ω.rec
    # split into left and right domains
    k = argmax(width(rec))
    Ψl = empty(Ψ)
    Ψr = empty(Ψ)
    for ψ in Ψ
        ψl, ψr = split(ψ, k)
        push!(Ψl, ψl)
        push!(Ψr, ψr)
    end
    Ω1 = MultiBernsteinCell(Ψl, copy(signs))
    Ω2 = MultiBernsteinCell(Ψr, copy(signs))
    return Ω1, Ω2
end

function restrict(Ω::MultiBernsteinCell{N,T}, k, surf) where {N,T}
    Ψ = Ω.Ψ
    ∇Ψ = Ω.∇Ψ
    signs = Ω.signs
    rec = Ω.rec
    xc = center(rec)
    Ψ̃ = BernsteinPolynomial{N - 1,T}[]
    ∇Ψ̃ = SVector{N - 1,BernsteinPolynomial{N - 1,T}}[] # one dimensional lower, so one less derivative in grad
    new_signs = empty(signs)
    for (ψ, s, ∇ψ) in zip(Ψ, signs, ∇Ψ)
        # why bound? dont we know that ∇Ψ[k] has a fixed sign on direction k?
        # pos_neg = bound(∇ψ[k],rec)[1] > 0 ? 1 : -1
        pos_neg = ∇ψ[k](xc) > 0 ? 1 : -1 # use sign?
        ψL = lower_restrict(ψ, k)
        sL = sgn(pos_neg, s, surf, -1)
        ψU = upper_restrict(ψ, k)
        sU = sgn(pos_neg, s, surf, 1)
        append!(Ψ̃, (ψL, ψU))
        append!(new_signs, (sL, sU))
    end
    Ω̃ = MultiBernsteinCell(Ψ̃, new_signs)
    return Ω̃
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

"""
    sgn(m,s,surface,side)

Compute the sign of the upper and lower restrictions of a function `ψ` with sign
`s`. The `side` variable is eithe `1` for the upper restriction, or `-1` for the
lower restriction, and `m` is the sign of `∇ψ ⋅ eₖ`.
"""
function sgn(m, s, surface, side)
    if m * s * side > 0 || surface
        if m * side > 0
            return 1
        else
            return -1
        end
    else
        return 0
    end
end

function StaticArrays.insert(x::Real, k, y::Real)
    return insert(SVector(x), k, y)
end
