#=
    Implementation of methods required to generate conforming elements from a
    level-set by recursive subdivision and one-dimensional root-finding.
=#

function quadgen(ls::LevelSet; qorder=5, maxdepth=20, maxslope=10, meshsize, porder=1)
    quadgen(CartesianLevelSet(ls; meshsize, order=porder);qorder, maxdepth, maxslope)
end

function quadgen(ls::CartesianLevelSet; qorder=5, maxdepth=20, maxslope=10)
    N = ambient_dimension(ls)
    qnodes = Vector{QuadratureNode{N,Float64}}()
    p   = (;qorder,maxdepth,maxslope)
    return quadgen!(qnodes, ls, p)
end

function quadgen!(qnodes, ls::CartesianLevelSet, p::NamedTuple)
    x1d, w1d = gausslegendre(p.qorder)
    s = levelset_sign(ls)
    # loop over the C∞ interpolants of `ls`
    for f in bernstein_interpolants(ls)
        root = MultiBernsteinCell([f], [s])
        surf = s == 0
        level = 0
        X, W  = _quadgen(root, surf, x1d, w1d, level, p)
        for (x, w) in zip(X, W) 
            n = ForwardDiff.gradient(f, x)
            n /= norm(n)
            κ = curvature(f, x)
            push!(qnodes, QuadratureNode(x, w, n, κ))
        end
    end        
    return qnodes
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

function dim1quad(Ψ::Vector{<:Function}, signs::Vector{<:Integer}, L, U, x1d, w1d, multi_zeros=true)
    roots = [L, U]
    if multi_zeros
        for ψ in Ψ
            union!(roots, find_zeros(ψ, L, U))
        end
    else
        for ψ in Ψ
            ψ(L) * ψ(U) < 0 && union!(roots, find_zero(ψ, (L, U)))
        end
    end
    sort!(roots)
    nodes   = Vector{Float64}()
    weights = Vector{Float64}()
    for (l, r) in zip(roots[1:end-1], roots[2:end])
        val = [ψ((l + r) / 2.0) for ψ in Ψ]
        if all((val .* signs) .≥ 0)
            # push then rescale the nodes
            append!(nodes, x1d)
            append!(weights, w1d)
            is = length(nodes)-length(x1d)+1 # where the new nodes start
            for i in is:length(nodes)
                nodes[i]   = nodes[i] * (r - l) + l
                weights[i] = weights[i] * (r - l)
            end
        end
    end
    return nodes, weights
end

function tensorquad(rec::HyperRectangle{D}, x1d, w1d) where {D}
    xl, xu = low_corner(rec), high_corner(rec)
    μ      = prod(xu-xl)
    nodes = map(Iterators.product(ntuple(i->x1d,D)...)) do x̂
        # map reference nodes to rec
        rec(SVector(x̂))
    end |> vec
    weights = map(Iterators.product(ntuple(i->w1d,D)...)) do ŵ
        # scale reference weights
        prod(ŵ)*μ
    end |> vec
    nodes, weights
end

function _quadgen(Ω::MultiBernsteinCell{N,T},surf,x1d, w1d, level, par) where {N,T}
    rec = Ω.rec
    xc  = center(rec)
    Ψ   = Ω.Ψ
    ∇Ψ  = Ω.∇Ψ
    signs  = Ω.signs
    D   = ambient_dimension(rec)
    @assert !surf || (D > 1)
    if level ≥ par.maxdepth
        @warn "Maximum depth reached: resorting to low-order quadrature"
        if surf
            return [xc], [prod(i->high_corner(rec)[i]-low_corner(rec)[i], D-1)]
        else
            return [xc], [prod(high_corner(rec).-low_corner(rec))]
        end
    end

    ctype = cell_type(Ω)
    # check for limiting cases of empty or whole cells
    if ctype == empty_cell
        return Vector{SVector{D,T}}(), Vector{T}()
    elseif ctype == whole_cell
        if surf
            return SVector{D,T}[], T[]
        else
            return tensorquad(rec,x1d,w1d)
        end
    end

    # base case
    if D == 1
        return dim1quad(Ψ, signs, low_corner(rec)[1], high_corner(rec)[1], x1d, w1d, true)
    end

    # find a heigh direction such that all of ∇Ψ are (provably) bounded away
    # from zero.
    bnds    = map(∇Ψ) do ∇ψ
        ntuple(d->bound(∇ψ[d],rec),D)
    end
    isvalid = ntuple(D) do dim
        all(bnds) do bnd
            (prod(bnd[dim])>0) &&
            (sum(bd->maximum(abs,bd), bnd)/minimum(abs,bnd[dim]) < par.maxslope)
        end
    end
    if !any(isvalid) # no valid direction so split
        Ω1,Ω2 = split(Ω)
        X1, W1 = _quadgen(Ω1, surf, x1d, w1d, level+1, par)
        X2, W2 = _quadgen(Ω2, surf, x1d, w1d, level+1, par)
        return (append!(X1, X2), append!(W1, W2))
    end

    # If there is a valid direction, we go down on it. Choose the direction which
    # is the least steep overall by maximizing the minimum of the derivative on
    # direction k over all functions
    ∇Ψc = map(∇Ψ) do ∇ψ
        ntuple(D) do d
            ∇ψc = abs.(∇ψ[d](xc))
            ∇ψc/norm(∇ψc)
        end
    end
    k = argmax(1:D) do dim
        if isvalid[dim]
            minimum(∇ψc -> abs(∇ψc[dim]),∇Ψc)
        else
            -Inf
        end
    end
    Ω̃  = restrict(Ω,k,surf) # the D-1 dimensional domain
    X, W = _quadgen(Ω̃, false, x1d, w1d, level, par)
    nodes = Vector{SVector{D,T}}()
    weights = Vector{T}()
    for (x, w) in zip(X, W)
        if surf
            # if we get here there should be only one Ψ
            @assert length(Ψ) == 1
            ψ = first(Ψ)
            ∇ψ = first(∇Ψ)
            lk, rk = low_corner(rec)[k], high_corner(rec)[k]
            ψₖ(y) = ψ(insert(x, k, y))
            if ψₖ(lk) * ψₖ(rk) < 0
                y = find_zero(ψₖ, (lk, rk))
                x̃ = insert(x, k, y)
                ∇ϕ = map(f->f(x̃),∇ψ)
                push!(nodes, x̃)
                push!(weights, w * norm(∇ϕ) / abs(∇ϕ[k]))
            end
        else
            Φ = [y -> ψ(insert(x, k, y)) for ψ in Ψ]
            Y, Ω = dim1quad(Φ, signs, low_corner(rec)[k], high_corner(rec)[k], x1d, w1d, false)
            for (y, ω) in zip(Y, Ω)
                push!(nodes, insert(x, k, y))
                push!(weights, w * ω)
            end
        end
    end
    ###########################################
    return nodes, weights
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
