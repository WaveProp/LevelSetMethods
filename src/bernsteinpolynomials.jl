"""
    struct BernsteinPolynomial{D,T}

Multivariate [Bernstein
polynomial](https://en.wikipedia.org/wiki/Bernstein_polynomial) given by tensor
product of one-dimesional Bernsteain polynomials on `domain`.

The exact expression is given by
```math
p(x_1,\\dots,x_D)=\\sum_{i_j=0}^{d_j}c_{i_1\\dots i_D}\\prod_{j=1}^D\\binom{d_j}{i_j}(x_j-l_j)^{i_j}(r_j-x_j)^{d_j-i_j}
```
where ``c_{i_1\\dots i_D}=\\texttt{coeffs}[i_1+1,\\dots,i_D+1]``, ``d_j=\\texttt{degree}[j]``,
``l_j=\\texttt{low\\_corner(domain)}[j]``,
``r_j=\\texttt{high\\_corner(domain)}[j]``
"""
struct BernsteinPolynomial{D,T} <: Function
    coeffs::Array{T,D}
    domain::HyperRectangle{D,T}
end

coeffs(p::BernsteinPolynomial) = p.coeffs
domain(p::BernsteinPolynomial) = p.domain
degree(p::BernsteinPolynomial) = size(coeffs(p)) .- 1

"""
    BernsteinPolynomial(c::Array,lc,uc)

Construct a [`BernsteinPolynomial`][@ref] with coefficients `c` on a
`HyperRectangle` with low corner `lc` and upper corner `uc`.
"""
function BernsteinPolynomial(c::Array{<:Any,D}, lc=ntuple(i -> 0.0, D),
                             uc=ntuple(i -> 1.0, D)) where {D}
    rec = HyperRectangle(lc, uc)
    return BernsteinPolynomial(c, rec)
end

"""
    lower_restrict(p::BernsteinPolynomial{D}, d::Integer) where{D}

Fast implementation of `partial_application(p, d, x)` when `x =
low_corner(p.domain)[d]`
"""
function lower_restrict(p::BernsteinPolynomial{D}, d::Integer) where {D}
    @assert 1 ≤ d ≤ D
    return BernsteinPolynomial(selectdim(p.coeffs, d, 1) .* 1, section(p.domain, d))
end

"""
    upper_restrict(p::BernsteinPolynomial{D}, d::Integer) where{D}

Fast implementation of `partial_application(p, d, x)` when `x = high_corner(p.domain)[d]`
"""
function upper_restrict(p::BernsteinPolynomial{D}, d::Integer) where {D}
    @assert 1 ≤ d ≤ D
    return BernsteinPolynomial(selectdim(p.coeffs, d, size(p.coeffs)[d]) .* 1,
                               section(p.domain, d))
end

@doc raw"""
    partial_application(p::BernsteinPolynomial{D,T},d::Integer, x::Real)

Return the `D-1`-variable Bernstein polynomial `p₁` by letting the `d`-th varible of `p` equal `x`, i.e.

```math
p_1(x_1,\dots,x_{d-1},x_{d+1},\dots,x_D)=p(x_1,\dots,x_{d-1},x,x_{d+1},\dots,x_D)
```
"""
function partial_application(p::BernsteinPolynomial{D,T}, d::Integer, x::Real) where {D,T}
    @assert 1 ≤ d ≤ D
    l = low_corner(p.domain)[d]
    r = high_corner(p.domain)[d]
    x = (x - l) / (r - l)
    sz = size(p.coeffs)
    sz′ = ntuple(i -> i < d ? sz[i] : sz[i + 1], D - 1)
    len = sz[d]
    out = similar(p.coeffs, sz′)
    for I in CartesianIndices(out)
        idxs = ntuple(i -> i < d ? (I[i]:I[i]) : i == d ? (1:len) : (I[i - 1]:I[i - 1]), D)
        c̃ = view(p.coeffs, idxs...)
        out[I] = evaluate_bernstein(SVector(x), c̃, Val{1}(), 1, len)
    end
    return BernsteinPolynomial(out,
                               ntuple(i -> i < d ? degree(p)[i] : degree(p)[i + 1], D - 1),
                               section(p.domain, d))
end

function (p::BernsteinPolynomial{D})(x::SVector{D}) where {D}
    l = low_corner(p.domain)
    r = high_corner(p.domain)
    x₀ = (x - l) ./ (r - l)
    return evaluate_bernstein(x₀, p.coeffs, Val{D}(), 1, length(p.coeffs))
end
(p::BernsteinPolynomial)(x) = p(SVector(x))

(P::SVector{N,<:BernsteinPolynomial})(x) where {N} = svector(i -> P[i](SVector(x)), N)

"""
    partials(p::BernsteinPolynomial)

Return a static vector containing the partial derivatives of `p` as
`BernsteinPolynomial`s.
"""
function partials(p::BernsteinPolynomial{D}, args...) where {D}
    svector(D) do d
        n = size(p.coeffs)[d]
        k = degree(p)[d]
        l = low_corner(p.domain)[d]
        u = high_corner(p.domain)[d]
        coeffs = mapslices(p.coeffs; dims=d) do b
            c = (b[2:n] .- b[1:(n - 1)]) .* k ./ (u - l)
            if k ≥ n
                push!(c, -b[n] * k / (u - l))
            end
            return c
        end
        return BernsteinPolynomial(coeffs, p.domain)
    end
end

function Base.show(io::IO, ::MIME"text/plain", p::BernsteinPolynomial)
    lb = low_corner(p.domain)
    ub = high_corner(p.domain)
    print(io, "Bernestein polynomial of degree ", degree(p), " on ",
          '[', lb[1], ',', ub[1], ']')
    for i in 2:length(lb)
        print(io, " × [", lb[i], ',', ub[i], ']')
    end
end

"""
    bound(p::BernsteinPolynomial)

Return the (approximated) lower and upper bound of `p` over `p.domain`
"""
function bound(p::BernsteinPolynomial, U::HyperRectangle = domain(p))
    @assert domain(p) == U
    extrema(coeffs(p))
end

"""
    Base.split(p::BernsteinPolynomial{D,T}, d::Integer, α=0.5)

Split `p.domain` along the `d`-th axis by a portion `0 ≤ α ≤ 1`.
Return the Bernstein polynomials on the two resulting domains.
"""
function Base.split(p::BernsteinPolynomial{D,T}, d::Integer, α=0.5) where {D,T}
    @assert 1 ≤ d ≤ D
    k = degree(p)[d]
    k == 0 && return p, p
    n = size(p.coeffs)[d]
    # FIXME: the code below is type-unstable. This is very likely due to
    # mapslices, which has some issues. Rewrite this as a standard loop to avoid
    # the use of high-order functions.
    coeffs = mapslices(p.coeffs; dims=d) do b
        c1 = Vector{T}()
        c2 = Vector{T}()
        for i in k:-1:1
            if i ≥ n
                push!(c1, b[1])
                @. b[1:(n - 1)] = b[1:(n - 1)] * (1 - α) + b[2:n] * α
                b[n] = b[n] * (1 - α)
            else
                push!(c1, b[1])
                pushfirst!(c2, b[i + 1])
                @. b[1:i] = b[1:i] * (1 - α) + b[2:(i + 1)] * α
            end
        end
        push!(c1, b[1])
        append!(c1, c2)
        return c1
    end
    split_point = low_corner(p.domain)[d] +
                  (high_corner(p.domain)[d] - low_corner(p.domain)[d]) * α
    rec1, rec2 = split(p.domain, d, split_point)
    p1 = BernsteinPolynomial(collect(selectdim(coeffs, d, 1:(k + 1))), rec1)
    p2 = BernsteinPolynomial(collect(selectdim(coeffs, d, (k + 1):(k + n))), rec2)
    return p1, p2
end

@doc raw"""
    rebase(a::Vector{<:Real}, l::Real, r::Real)

Given the vector of coefficients `a` of a polynomial in monomial basis, 
return the vector of coefficients `ã` in the basis after an affine transform.

```math
\sum_{i=0}^{n-1}a[i+1]x^i = \sum_{i=0}^{n-1}\tilde{a}[i+1](\frac{x-l}{r-l})^i
```
"""
function rebase(a::Vector{<:Real}, l::Real, r::Real)
    n = length(a)
    ã = copy(a)
    for i in 0:(n - 2)
        ã[(n - i):n] .*= (r - l)
        ã[(n - i - 1):(n - 1)] .+= ã[(n - i):n] ./ (r - l) .* l
    end
    return ã
end

function rebase(A::Array{<:Real,D}, rec::HyperRectangle{D}) where {(D)}
    L = low_corner(rec)
    R = high_corner(rec)
    Ã = copy(A)
    for d in 1:D
        Ã = mapslices(Ã; dims=d) do a
            return rebase(a, L[d], R[d])
        end
    end
    return Ã
end

"""
    power2bernstein(a::Array{<:Real,D}, U::HyperRectangle{D}) where{D}

Convert a polynomial in power series into a Bernstein polynomial on `U` of minimal degree.
"""
function power2bernstein(a::Array{<:Real,D}, U::HyperRectangle{D}) where {D}
    k = size(a) .- 1
    b = zeros(k .+ 1)
    a = rebase(a, U)
    for i in CartesianIndices(a)
        temp = zeros(Tuple([i[j] for j in 1:D]))
        for l in CartesianIndices(temp)
            temp[l] = a[l] * prod(1:D) do j
                             return binomial(i[j] - 1, l[j] - 1) / binomial(k[j], l[j] - 1)
                             end
        end
        b[i] = sum(temp)
    end
    return BernsteinPolynomial(b, U)
end

"""
    raise_degree(p::BernsteinPolynomial{D,T}, k::NTuple{D,Int}) where {D,T}

Raise the degree of the BernsteinPolynomial `p` to `k`.
"""
function raise_degree(p::BernsteinPolynomial{D,T}, k::NTuple{D,Int}) where {D,T}
    if !all(k .>= degree(p))
        @error "New degree should be higher."
        return p
    end
    C = zeros(T, k.+1)
    n = size(coeffs(p))
    inds = CartesianIndices(ntuple(d->1:n[d],D))
    copyto!(C, inds, coeffs(p), inds)
    for d in 1:D 
        for i in 1:k[d]-n[d]+1
            selectdim(C, d, 1+i:n[d]+i) .+= selectdim(C, d, i:n[d]+i-1)
        end
    end
    return BernsteinPolynomial(C, domain(p))
end


# Tensor-product pattern adapted from FastChebInterp.jl (MIT license)
@fastmath function evaluate_bernstein(x::SVector{N}, c::AbstractArray, ::Val{dim}, i1,
                            len) where {N,dim}
    n = size(c, dim)
    @inbounds xd = x[dim]
    # see https://personal.math.ubc.ca/~cass/graphics/text/www/pdf/a6.pdf for
    # one dimensional case
    if dim == 1
        s = 1 - xd
        @inbounds P = c[i1]
        C = (n - 1) * xd
        for k in 1:(n - 1)
            @inbounds P = P * s + C * c[i1 + k]
            C = C * (n - k - 1) / (k + 1) * xd
        end
        return P
    else
        Δi = len ÷ n # column-major stride of current dimension

        # we recurse downward on dim for cache locality,
        # since earlier dimensions are contiguous
        dim′ = Val{dim - 1}()

        s = 1 - xd
        P = evaluate_bernstein(x, c, dim′, i1, Δi)
        C = (n - 1) * xd
        for k in 1:(n - 1)
            P = P * s + C * evaluate_bernstein(x, c, dim′, i1 + k * Δi, Δi)
            C = C * (n - k - 1) / (k + 1) * xd
        end
        return P
    end
end

# deCasteljau algorithm for splitting a Bernstein polynomial
function deCasteljau!(coeffs, d, t)
    N = ndims(coeffs)
    sz = size(coeffs)
    L   = LinearIndices(sz)
    I   = ntuple(i->i==d ? (1:1) : 1:sz[i], N)
    pivots = view(L,I...)
    Δl  = 1
    for i in 1:d-1
        Δl *= sz[i]
    end
    @assert 0 ≤ t ≤ 1
    n   = size(coeffs,d) - 1
    for shift in pivots
        for i in 1:n
            val = coeffs[shift]
            jend= size(coeffs,d) - i
            for j in 1:size(coeffs,d)-i
                @inbounds coeffs[shift+(j-1)*Δl] = (1-t)*coeffs[shift+(j-1)*Δl] + t*coeffs[shift+j*Δl]
            end
            coeffs[shift+jend*Δl] = val
        end
    end
    return coeffs
end
deCasteljau(coeffs,d,t) = deCasteljau!(copy(coeffs),d,t)



@enum CellType empty_cell whole_cell cut_cell
# Level set expressed by a Vector of BernsteinPolynomial
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

@doc raw"""
    curvature(p::BernsteinPolynomial{D}, x::SVector{D})

Calculate the curvature of the level-set of `p` at `x` by the following formula

```math
\kappa = \nabla\cdot\frac{\nabla p}{|\nabla p|} = \frac{\Delta p}{|\nabla p|} + \frac{\nabla^\perp p\nabla^2 p\nabla p}{|\nabla p|^3}
```
"""
function curvature(p::BernsteinPolynomial{D}, x::SVector{D}) where {D}
    ∇p = ForwardDiff.gradient(p, x)
    Hp = ForwardDiff.hessian(p, x)
    Np = norm(∇p)

    tr(Hp)/Np - dot(∇p, Hp, ∇p)/Np^3
end
