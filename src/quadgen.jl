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
    x1d .= (x1d .+ 1) ./ 2
    w1d .= w1d ./ 2
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