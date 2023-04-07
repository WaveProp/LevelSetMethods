"""
    bound(f,l,u)

Return a tuple `(lb,ub)` of a rigorous (but not necessarily sharp) lower and
upper bound for the function `f` on the axis-aligned hyperrectangle with
low-corner `l` and upper corner `u`.

By default, uses `IntervalArithmetic` to compute the bounds. Overload this
function for your own types if a more efficient/accurate `bound` function is
known.

```jldoctest
bound(x->x[1],-1,1)

# output
(-1.0, 1.0)
```

```jldoctest
bound(x->x[1]^2,-1,1)

# output
(0.0, 1.0)
```

"""
function bound(f,l::SVector{N,<:Real},u::SVector{N,<:Real}) where {N}
    @assert length(l) === length(u) "l and u must have the same length"
    box = IntervalBox(ntuple(i->l[i]..u[i],N))
    I = f(box)
    return I.lo, I.hi
end
bound(f,l,u) = bound(f,SVector(l),SVector(u))
bound(f,rec::WPB.HyperRectangle) = bound(f, WPB.low_corner(rec), WPB.high_corner(rec))


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

function my_find_zero(f,(l,u))
    if min(abs(f(l)),abs(f(u))) < 1e-8
        return abs(f(l)) < abs(f(u)) ? l : u
    elseif f(l)*f(u) < 0
        return find_zero(f,(l,u))
    else
        return nothing
    end
end

function neighbors(I::CartesianIndex{D}, size::NTuple{D}, occupied=Set{CartesianIndex{D}}()) where {D}
    N = Vector{CartesianIndex{D}}()
    for d in 1:D
        if I[d] > 1
            J = CartesianIndex(ntuple(i->i == d ? I[i]-1 : I[i], D))
            !(J in occupied) && push!(N, J)
        end
        if I[d] < size[d]
            J = CartesianIndex(ntuple(i->i == d ? I[i]+1 : I[i], D))
            !(J in occupied) && push!(N, J)
        end
    end
    return N
end