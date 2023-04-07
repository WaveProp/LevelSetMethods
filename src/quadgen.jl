#=
    Implementation of methods required to generate quadratures from a
    level-set by recursive subdivision and one-dimensional root-finding.
=#

function quadgen(ls::LevelSet; qorder=5, maxdepth=20, maxslope=10, curvature=false, merge=false)
    N = ambient_dimension(ls)
    M = geometric_dimension(ls)
    qrule1d = WPB.qrule_for_reference_shape(WPB.ReferenceLine(), qorder)
    nq = length(qrule1d()[2])^M
    qnodes = Vector{QuadratureNode{N,Float64}}()
    p   = (;qrule1d,maxdepth,maxslope,qorder,curvature,merge)
    qnodes, F = quadgen!(qnodes, ls, p)
    return reshape(qnodes, nq, :), F
end

function quadgen!(qnodes, ls::LevelSet, p::NamedTuple)
    x1d,w1d = p.qrule1d()
    x1d = [x[1] for x in x1d] |> Vector
    w1d = [w[1] for w in w1d] |> Vector
    s = levelset_sign(ls)
    ϕ = levelset_function(ls)

    if p.merge
        # first loop: estimate the intersection volume in each Cartesian cell
        interps = interpolants(ϕ)
        volume = Dict{CartesianIndex, Float64}()
        surf = s == 0
        _p = (maxdepth=5,maxslope=10,qorder=3)
        for (I, f) in interps
            root = MultiBernsteinCell([f], [s])
            _, W = _quadgen(root, surf, x1d, w1d, 0, _p)
            volume[I] = sum(W)
        end

        # second loop: find the cells with small intersection volume and merge with biggest neighbor
        F = Vector{BernsteinPolynomial}()
        msh = mesh(ϕ); stp = step(msh); sz = size(msh); D = ambient_dimension(ϕ)
        small_cell   = Set{CartesianIndex}()
        if surf
            threshold = prod(stp) / maximum(stp) / 2^(D+1)
        end
        for (I, v) in volume
            0 < v < threshold && push!(small_cell, I)
        end
        treated_cell = Set{CartesianIndex}()
        if surf
            for I in small_cell
                neighs = neighbors(I, sz, treated_cell)
                neigh_vol = [volume[J] for J in neighs]
                if isempty(neigh_vol)
                    @warn "Some small cells can't be merged"
                    push!(F, interpolant(ϕ, I))
                    continue
                end
                J = neighs[argmax(neigh_vol)]
                push!(treated_cell, I, J)
                push!(F, merge_interpolant(ϕ, I, J))
            end
            for (I, _) in volume
                if !(I in treated_cell || I in small_cell)
                    push!(treated_cell, I)
                    push!(F, interpolant(ϕ, I))
                end
            end
        end
    else
        F = [interp[2] for interp in interpolants(ϕ)]
    end

    # loop over the C∞ interpolants of `ls`
    for f in F
        root = MultiBernsteinCell([f], [s])
        surf = s == 0
        level = 0
        X, W  = _quadgen(root, surf, x1d, w1d, level, p)
        for (x, w) in zip(X, W)
            n = ForwardDiff.gradient(f, x)
            n /= norm(n)
            κ = p.curvature ? curvature(f, x) : nothing
            push!(qnodes, QuadratureNode(x, w, n, κ))
        end
    end
    return qnodes, F
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
            # return [xc], [prod(i->high_corner(rec)[i]-low_corner(rec)[i], D-1)]
            return [], []
        else
            # return [xc], [prod(high_corner(rec).-low_corner(rec))]
            return [], []
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
        nodes, weights = dim1quad(Ψ, signs, low_corner(rec)[1], high_corner(rec)[1], x1d, w1d, true)
        return nodes,weights
    end

    # find candidate heigh directions such that all of ∇Ψ are (provably) bounded away
    # from zero.
    candidates = []
    partial_bnds = [[bound(∇ψ[d],rec) for d in 1:D] for ∇ψ in ∇Ψ]
    grad_l1_norm = [sum(x->maximum(abs.(x)), row) for row in partial_bnds]
    for d in 1:D
        cand = true
        for j in 1:length(∇Ψ)
            if prod(partial_bnds[j][d]) ≤ 0 || grad_l1_norm[j] ≥ minimum(abs.(partial_bnds[j][d]))*par.maxslope
                cand = false; break
            end
        end
        if cand
            push!(candidates, d)
        end
    end
    if isempty(candidates) # no valid direction so split
        Ω1,Ω2 = split(Ω)
        X1, W1 = _quadgen(Ω1, surf, x1d, w1d, level+1, par)
        X2, W2 = _quadgen(Ω2, surf, x1d, w1d, level+1, par)
        return (append!(X1, X2), append!(W1, W2))
    end

    # If there is a valid direction, we go down on it. Choose the direction with
    # the least steep overall by maximizing the minimum of the derivative on
    # direction k over all functions
    partial_midvals  = [[abs(∇ψ[i](xc)) for i in 1:D] for ∇ψ in ∇Ψ]
    grad_l1_norm_mid = [sum(row) for row in partial_midvals]
    score = [sum(j->partial_midvals[j][i]/grad_l1_norm_mid[j], 1:length(∇Ψ)) for i in candidates]
    k = candidates[argmax(score)]

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
            y = my_find_zero(ψₖ, (lk, rk))
            isnothing(y) && continue
            x̃ = insert(x, k, y)
            ∇ϕ = map(f->f(x̃),∇ψ)
            push!(nodes, x̃)
            push!(weights, w * norm(∇ϕ) / abs(∇ϕ[k]))
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

function dim1quad(Ψ::Vector{<:Function}, signs::Vector{<:Integer}, L, U, x1d, w1d, multi_zeros=true;tol=10^(-12))
    roots = [L, U]
    if multi_zeros
        for ψ in Ψ
            for z in find_zeros(ψ, L, U)
                if all(r->abs(r-z)>tol, roots)
                    push!(roots, z)
                end
            end
        end
    else
        for ψ in Ψ
            y = my_find_zero(ψ, (L, U))
            isnothing(y) || union!(roots, y)
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

function NystromMesh(ls::LevelSet;kwargs...)
    N = ambient_dimension(ls)
    qnodes, F = quadgen(ls;kwargs...)
    nq,nel = size(qnodes)
    etype2qtags = Dict{DataType,Matrix{Int}}(Nothing => reshape(1:length(qnodes),nq,:))
    return NystromMesh{N,Float64}(;qnodes=vec(qnodes),etype2qtags)
end