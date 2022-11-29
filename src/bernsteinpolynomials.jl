"""
    struct Polynomial{D,T}

`D`-dimensional polynomial with coefficients of type `T`. The coefficients are
stored as a dense array, and are implicitly associated with a monomial basis;
the `c[I]` coefficient multiplies the monomial term `prod(x.^(I .- 1))`, where
`I` is a `D`-dimensional multi-index.
"""
struct Polynomial{D,T}
    coeffs::Array{T,D}
end

coeffs(p::Polynomial) = p.coeffs

const VARIABLE_LABELS = ["x₁" "x₂" "x₃"]
const POWER_LABELS    = ["" "²" "³" "⁴" "⁵" "⁶"]

variable_label(i::Int) = i > length(VARIABLE_LABELS) ? "x[$i]" : VARIABLE_LABELS[i]
power_label(i::Int)    = i > length(POWER_LABELS)    ? "^$i" : POWER_LABELS[i]

function Base.show(io::IO, ::MIME"text/plain",p::Polynomial{D}) where {D}
    str = ""
    C   = p.coeffs
    for I in CartesianIndices(C)
        c = C[I]
        iszero(c) && continue
        θ = Tuple(I) .- 1
        m = [iszero(θ[d]) ? "" : "$(variable_label(d))$(power_label(θ[d]))" for d in 1:D]
        times_str = sum(Tuple(I)) == D ? "" : "⋅"
        str = str * "$c" * times_str * prod(m) * " + "
    end
    print(io, str[1:end-2])
end

# TODO: evalate `Polynomial` using a multidimensional Horner algorithm
function (p::Polynomial{N,T})(x::SVector{N,S}) where {N,T,S}
    c = coeffs(p)
    return sum(CartesianIndices(c)) do I
        k = Tuple(I) .- 1
        c[I] * prod(x .^ k)
    end
end

(p::Polynomial)(x) = p(SVector(x))

(P::SVector{N,<:Polynomial})(x) where {N} = svector(i->P[i](SVector(x)), N)

"""
    struct BernsteinPolynomial{D,T}
        coeffs::Array{T,D}
        degree::NTuple{D,Integer}
        domain::HyperRectangle{D,T}
    end

Create a multivariate Bernstein polynomial

```math
p(x_1,\\dots,x_D)=\\sum_{i_j=0}^{d_j}c_{i_1\\dots i_D}\\prod_{j=1}^D\\binom{d_j}{i_j}(x_j-l_j)^{i_j}(r_j-x_j)^{d_j-i_j}
```
where ``c_{i_1\\dots i_D}=\\texttt{coeffs}[i_1+1,\\dots,i_D+1]``, ``d_j=\\texttt{degree}[j]``,
``l_j=\\texttt{low\\_corner(domain)}[j]``, ``r_j=\\texttt{high\\_corner(domain)}[j]``

"""
struct BernsteinPolynomial{D,T} <: Function
    coeffs::Array{T,D}
    degree::NTuple{D,Integer}
    domain::HyperRectangle{D,T}
end

order(p::BernsteinPolynomial) = p.degree
domain(p::BernsteinPolynomial) = p.domain

"""
    reference_cube(::Val{D})

Return the `[0,1]ᴰ` reference domain as a `HyperRectangle{D,Float64}`.
"""
reference_cube(::Val{D}) where {D} = HyperRectangle(svector(i->0., D), svector(i->1., D))
reference_cube(D) = HyperRectangle(svector(i->0., D), svector(i->1., D))

BernsteinPolynomial(c::Array{<:Real,D}) where {D} = BernsteinPolynomial(c, size(c).-1, reference_cube(Val(D)))


"""
    lower_restrict(p::BernsteinPolynomial{D}, d::Integer) where{D}

Fast implemtation of `partial_application(p, d, x)` when `x = low_corner(p.domain)[d]`
"""
function lower_restrict(p::BernsteinPolynomial{D}, d::Integer) where{D}
    @assert 1 ≤ d ≤ D
    BernsteinPolynomial(selectdim(p.coeffs, d, 1).*1, ntuple(i->i<d ? p.degree[i] : p.degree[i+1], D-1), section(p.domain, d))
end

"""
    upper_restrict(p::BernsteinPolynomial{D}, d::Integer) where{D}

Fast implemtation of `partial_application(p, d, x)` when `x = high_corner(p.domain)[d]`
"""
function upper_restrict(p::BernsteinPolynomial{D}, d::Integer) where{D}
    @assert 1 ≤ d ≤ D
    if p.degree[d] ≥ size(p.coeffs)[d]
        BernsteinPolynomial(reshape([0.], ntuple(i->1,D-1)), ntuple(i->1,D-1), section(p.domain, d))
    else
        BernsteinPolynomial(selectdim(p.coeffs, d, size(p.coeffs)[d]).*1, ntuple(i->i<d ? p.degree[i] : p.degree[i+1], D-1), section(p.domain, d))
    end
end

@doc raw"""
    partial_application(p::BernsteinPolynomial{D,T},d::Integer, x::Real)

Return the `D-1`-variable Bernstein polynomial `p₁` by letting the `d`-th varible of `p` equal `x`, i.e.

```math
p_1(x_1,\dots,x_{d-1},x_{d+1},\dots,x_D)=p(x_1,\dots,x_{d-1},x,x_{d+1},\dots,x_D)
```
"""
function partial_application(p::BernsteinPolynomial{D,T},d::Integer, x::Real) where{D,T}
    @assert 1 ≤ d ≤ D
    l = low_corner(p.domain)[d]; r = high_corner(p.domain)[d]
    x = (x - l) / (r - l)
    sz = size(p.coeffs)
    sz′ =  ntuple(i -> i < d ? sz[i] : sz[i+1], D-1)
    len =  sz[d]
    out = similar(p.coeffs,sz′)
    for I in CartesianIndices(out)
        idxs = ntuple(i->i < d ? (I[i]:I[i]) : i==d ? (1:len) : (I[i-1]:I[i-1]),D)
        c̃    = view(p.coeffs,idxs...)
        out[I] = evaluate(SVector(x),c̃,Val{1}(),1,len)
    end
    BernsteinPolynomial(out, ntuple(i->i<d ? p.degree[i] : p.degree[i+1], D-1), section(p.domain, d))
end

function (p::BernsteinPolynomial{D})(x::SVector{D}) where{D}
    l = low_corner(p.domain); r = high_corner(p.domain)
    x₀ = (x - l) ./ (r - l)
    evaluate(x₀, p.coeffs, Val{D}(), 1, length(p.coeffs))
end
(p::BernsteinPolynomial)(x) = p(SVector(x))

(P::SVector{N,<:BernsteinPolynomial})(x) where {N} = svector(i->P[i](SVector(x)), N)

"""
    gradient(p::BernsteinPolynomial{D}) where{D}

Return the gradien of `p` as an `SVector{D,BernsteinPolynomial{D}}`
"""
function gradient(p::BernsteinPolynomial{D}) where{D}
    svector(D) do d
        n = size(p.coeffs)[d]
        k = p.degree[d]
        l = low_corner(p.domain)[d]; u = high_corner(p.domain)[d]
        coeffs = mapslices(p.coeffs, dims=d) do b
            c = (b[2:n] .- b[1:n-1]) .* k ./ (u-l)
            if k ≥ n
                push!(c, -b[n]*k/(u-l))
            end
            c
        end
        BernsteinPolynomial(coeffs, ntuple(i->i==d ? p.degree[i]-1 : p.degree[i], D), p.domain)
    end
end

function Base.show(io::IO, ::MIME"text/plain", p::BernsteinPolynomial)
    lb = low_corner(p.domain)
    ub = high_corner(p.domain)
    print(io, "Bernestein polynomial of order ", order(p), " on ",
          '[', lb[1], ',', ub[1], ']')
    for i = 2:length(lb)
        print(io, " × [", lb[i], ',', ub[i], ']')
    end
end

"""
    bound(p::BernsteinPolynomial)

Return the (approximated) lower and upper bound of `p` over `p.domain`
"""
function bound(p::BernsteinPolynomial, rec::HyperRectangle=domain(p))
    # @assert domain(p) == rec # maybe it makes sense to rescale p if domain(p) is not rec?
    M = maximum(p.coeffs)
    m = minimum(p.coeffs)
    if M < 0 && prod(size(p.coeffs)) < prod(p.degree.+1)
        M = 0
    end
    if m > 0 && prod(size(p.coeffs)) < prod(p.degree.+1)
        m = 0
    end
    m, M
end

"""
    Base.split(p::BernsteinPolynomial{D,T}, d::Integer, α=0.5)

Split `p.domain` along the `d`-th axis by a portion `0 ≤ α ≤ 1`.
Return the Bernstein polynomials on the two resulting domains.
"""
function Base.split(p::BernsteinPolynomial{D,T}, d::Integer, α=0.5) where {D,T}
    @assert 1 ≤ d ≤ D
    k = p.degree[d]
    k == 0 && return p, p
    n = size(p.coeffs)[d]
    coeffs = mapslices(p.coeffs, dims=d) do b
        c1 = Vector{T}(); c2 = Vector{T}()
        for i in k:-1:1
            if i ≥ n
                push!(c1, b[1])
                @. b[1:n-1] = b[1:n-1]*(1-α) + b[2:n]*α
                b[n] = b[n]*(1-α)
            else
                push!(c1, b[1])
                pushfirst!(c2, b[i+1])
                @. b[1:i] = b[1:i]*(1-α) + b[2:i+1]*α
            end
        end
        push!(c1, b[1])
        append!(c1, c2)
        c1
    end
    split_point = low_corner(p.domain)[d] + (high_corner(p.domain)[d] - low_corner(p.domain)[d])*α
    rec1, rec2 = split(p.domain, d, split_point)
    p1 = BernsteinPolynomial(collect(selectdim(coeffs, d, 1:k+1)), p.degree, rec1)
    p2 = BernsteinPolynomial(collect(selectdim(coeffs, d, k+1:k+n)), p.degree, rec2)
    p1, p2
end

"""
    Cartesian_grid(p::BernsteinPolynomial{D,T}, grids::NTuple{D,Integer}) where {D,T}

Divide `p.domain` into an equidistant `grids[1]×⋯×grids[D]` Cartesian grid.
Return the list of Bernstein polynomials on each grid cell.
"""
function Cartesian_grid(p::BernsteinPolynomial{D,T}, grids::NTuple{D,Integer}) where {D,T}
    list_final = [p]
    list_tempo = Vector{BernsteinPolynomial{D,T}}()
    for d in D:-1:1
        if grids[d] > 1
            for q in list_final
                q2 = q
                for g in grids[d]:-1:2
                    q1, q2 = split(q2, d, 1/g)
                    push!(list_tempo, q1)
                end
                push!(list_tempo, q2)
            end
            list_final = list_tempo
            list_tempo = Vector{BernsteinPolynomial{D,T}}()
        end
    end
    list_final
end

"""
    Cartesian_grid(p::BernsteinPolynomial{D,T}, grid::Integer) where {D,T}

Simple interface for `Cartesian_grid(p, ntuple(i->grid, D))`
"""
Cartesian_grid(p::BernsteinPolynomial{D,T}, grid::Integer) where {D,T} = Cartesian_grid(p, ntuple(i->grid
, D))

function rebase(a::Vector{<:Real}, l::Real, r::Real)
    n = length(a)
    ã = copy(a)
    for i in 0:n-2
        ã[n-i:n] .*= (r-l)
        ã[n-i-1:n-1] .+= ã[n-i:n] ./ (r-l) .* l
    end
    ã
end

function rebase(A::Array{<:Real,D}, rec::HyperRectangle{D}) where(D)
    L = low_corner(rec)
    R = high_corner(rec)
    Ã = copy(A)
    for d in 1:D
        Ã = mapslices(Ã, dims=d) do a
            rebase(a, L[d], R[d])
        end
    end
    Ã
end

"""
    power2bernstein(a::Array{<:Real,D}, U::HyperRectangle{D}=□(D), k=size(a).-1) where{D}

Convert a polynomial in power series into a Bernstein polynomial on `U` of degree `k`.
"""
function power2bernstein(a::Array{<:Real,D}, U::HyperRectangle{D}=reference_cube(D), k=size(a).-1) where{D}
    b = zeros(k.+1)
    a = rebase(a, U)
    for i in CartesianIndices(a)
        temp = zeros(Tuple([i[j] for j = 1:D]))
        for l in CartesianIndices(temp)
            temp[l] = a[l] * prod(1:D) do j
                binomial(i[j]-1, l[j]-1) / binomial(k[j], l[j]-1)
            end
        end
        b[i] = sum(temp)
    end
    BernsteinPolynomial(b, k, U)
end

# Adapted from FastChebInterp
function evaluate(x::SVector{N}, c::AbstractArray, ::Val{dim}, i1, len) where {N,dim}
    n = size(c,dim)
    @inbounds xd = x[dim]
    # idea taken from here https://personal.math.ubc.ca/~cass/graphics/text/www/pdf/a6.pdf
    if dim == 1
        s = 1-xd
        @inbounds P = c[i1]
        C = (n-1)*xd
        for k in 1:n-1
            @inbounds P = P*s + C*c[i1+k]
            C = C*(n-k-1)/(k+1)*xd
        end
        return P
    else
        Δi = len ÷ n # column-major stride of current dimension

        # we recurse downward on dim for cache locality,
        # since earlier dimensions are contiguous
        dim′ = Val{dim-1}()

        s = 1-xd
        P = evaluate(x, c, dim′, i1, Δi)
        C = (n-1)*xd
        for k in 1:n-1
            P = P*s + C*evaluate(x,c,dim′,i1+k*Δi,Δi)
            C = C*(n-k-1)/(k+1)*xd
        end
        return P
    end
end
