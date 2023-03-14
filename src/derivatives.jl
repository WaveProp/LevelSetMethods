"""
    abstract type SpatialScheme

Discretization schemes for the spatial derivatives.
"""
abstract type SpatialScheme end

"""
    struct Upwind <: SpatialScheme

First-order upwind method. Uses `D⁺` and `D⁻` to compute the upwind derivative
approximation.
"""
struct Upwind <: SpatialScheme end

"""
    struct WENO5 <: SpatialScheme

Fith-order weighted essentially non-osciallatory method. Uses [`weno5⁺`](@ref) and
[`weno5⁻`](@ref) to compute right and left biased derivaties.
"""
struct WENO5 <: SpatialScheme end

"""
    D⁰(ϕ,I,dim)

Centered finite difference scheme for first order derivative at grid point `I`
along dimension `dim`.
"""
function D⁰(ϕ,I,dim,h=step(ϕ,dim))
    Im   = decrement_index(I,dim)
    Ip   = increment_index(I,dim)
    return (ϕ[Ip] - ϕ[Im]) / (2h)
end

"""
    D⁺(ϕ,I,dim)

Forward finite difference scheme for first order derivative at grid point `I`
along dimension `dim`.
"""
@inline function D⁺(ϕ,I,dim,h=step(ϕ,dim))
    Ip = increment_index(I,dim)
    return (ϕ[Ip] - ϕ[I]) / h
end

function D⁺⁺(ϕ,I,dim,h=step(ϕ,dim))
    Ip   = increment_index(I,dim)
    Ipp  = increment_index(I,dim,2)
    return (-1.5*ϕ[I] + 2*ϕ[Ip] - 1/2*ϕ[Ipp]) / h
end

"""
    D⁻(ϕ,I,dim)

Backward finite difference scheme for first order derivative at grid point `I`
along dimension `dim`.
"""
function D⁻(ϕ,I,dim,h=step(ϕ,dim))
    Im   = decrement_index(I,dim)
    return (ϕ[I] - ϕ[Im]) / h
end

function D⁻⁻(ϕ,I,dim,h=step(ϕ,dim))
    Im   = decrement_index(I,dim)
    Imm   = decrement_index(I,dim,2)
    return (1.5*ϕ[I] - 2*ϕ[Im] + 1/2*ϕ[Imm]) / h
end

"""
    weno5⁻(ϕ,I,dim[, h])

Fith-order weno derivative using the stencil `i-3,i-2,i-1,i,i+1,i+2`.
"""
function weno5⁻(ϕ,I,dim,h=step(ϕ,dim))
    # see section 3.4 of Osher-Fedwik
    Im  = decrement_index(I,dim)
    Imm = decrement_index(Im,dim)
    Ip  = increment_index(I,dim)
    Ipp = increment_index(Ip,dim)
    # finite differences
    v1  = D⁻(ϕ,Imm,dim,h)
    v2  = D⁻(ϕ,Im,dim,h)
    v3  = D⁻(ϕ,I,dim,h)
    v4  = D⁻(ϕ,Ip,dim,h)
    v5  = D⁻(ϕ,Ipp,dim,h)
    # third-order estimates
    dϕ1 =  (1/3)*v1 - (7/6)*v2 + (11/6)*v3
    dϕ2 = -(1/6)*v2 + (5/6)*v3 + (1/3)*v4
    dϕ3 = (1/3)*v3 + (5/6)*v4 - (1/6)*v5
    # smoothness indicators
    S1 = (13/12)*(v1 - 2*v2 + v3)^2 + (1/4)*(v1 - 4*v2 + 3*v3)^2
    S2 = (13/12)*(v2 - 2*v3 + v4)^2 + (1/4)*(v2 - v4)^2
    S3 = (13/12)*(v3 - 2*v4 + v5)^2 + (1/4)*(3*v3 - 4*v4 + v5^3)^2
    # fudge factor
    ϵ = 1e-6*max(v1^2,v2^2,v3^2,v4^2,v5^2) + 1e-99
    # weights
    α1 = 0.1 / (S1+ϵ)^2
    α2 = 0.6 / (S2+ϵ)^2
    α3 = 0.3 / (S3+ϵ)^2
    ω1 = α1 / (α1 + α2 + α3)
    ω2 = α2 / (α1 + α2 + α3)
    ω3 = α3 / (α1 + α2 + α3)
    # WENO approximation
    return ω1*dϕ1 + ω2*dϕ2 + ω3*dϕ3
end

"""
    weno5⁺(ϕ,I,dim[,h])

Fith-order weno derivative using the stencil `i-2,i-1,i,i+1,i+2,i+3`.
"""
function weno5⁺(ϕ,I,dim,h=step(ϕ,dim))
    # see section 3.4 of Osher-Fedwik
    Im  = decrement_index(I,dim)
    Imm = decrement_index(Im,dim)
    Ip  = increment_index(I,dim)
    Ipp = increment_index(Ip,dim)
    # finite differences
    v1  = D⁺(ϕ,Ipp,dim,h)
    v2  = D⁺(ϕ,Ip,dim,h)
    v3  = D⁺(ϕ,I,dim,h)
    v4  = D⁺(ϕ,Im,dim,h)
    v5  = D⁺(ϕ,Imm,dim,h)
    # third-order estimates
    dϕ1 =  (1/3)*v1 - (7/6)*v2 + (11/6)*v3
    dϕ2 = -(1/6)*v2 + (5/6)*v3 + (1/3)*v4
    dϕ3 = (1/3)*v3 + (5/6)*v4 - (1/6)*v5
    # smoothness indicators
    S1 = (13/12)*(v1 - 2*v2 + v3)^2 + (1/4)*(v1 - 4*v2 + 3*v3)^2
    S2 = (13/12)*(v2 - 2*v3 + v4)^2 + (1/4)*(v2 - v4)^2
    S3 = (13/12)*(v3 - 2*v4 + v5)^2 + (1/4)*(3*v3 - 4*v4 + v5^3)^2
    # fudge factor
    ϵ = 1e-6*max(v1^2,v2^2,v3^2,v4^2,v5^2) + 1e-99
    # weights
    α1 = 0.1 / (S1+ϵ)^2
    α2 = 0.6 / (S2+ϵ)^2
    α3 = 0.3 / (S3+ϵ)^2
    ω1 = α1 / (α1 + α2 + α3)
    ω2 = α2 / (α1 + α2 + α3)
    ω3 = α3 / (α1 + α2 + α3)
    # WENO approximation
    return ω1*dϕ1 + ω2*dϕ2 + ω3*dϕ3
end

"""
    D2⁰(ϕ,I,dim[, h])

Centered finite difference scheme for second order derivative at grid point `I`
along dimension `dim`. E.g. if `dim=1`, this approximates `∂ₓₓ`.
"""
function D2⁰(ϕ,I,dim,h=step(ϕ,dim))
    Im   = decrement_index(I,dim)
    Ip   = increment_index(I,dim)
    return (ϕ[Ip] - 2ϕ[I] + ϕ[Im]) / h^2
end

"""
    D2(ϕ,I,dims)

Finite difference scheme for second order derivative at grid point `I`
along the dimensions `dims`.

If `dims[1] == dims[2]`, it is more efficient to call `D2⁰(ϕ,I,dims[1])`.
"""
function D2(ϕ,I,dims)
    h = map(d->step(ϕ,d),dims)
    Ip = increment_index(I,dims[1])
    Im = decrement_index(I,dims[1])
    (D⁰(ϕ,Ip,dims[2]) - D⁰(ϕ,Im,dims[2]))/(2*h[dims[1]])
end

"""
    D2⁺⁺(ϕ,I,dim)

Upward finite difference scheme for second order derivative at grid point `I`
along dimension `dim`. E.g. if `dim=1`, this approximates `∂ₓₓ`.
"""
function D2⁺⁺(ϕ,I,dim,h=step(ϕ,dim))
    Ip   = increment_index(I,dim,1)
    Ipp  = increment_index(I,dim,2)
    return (ϕ[I] - 2ϕ[Ip] + ϕ[Ipp]) / h^2
end

"""
    D2⁻⁻(ϕ,I,dim)

Backward finite difference scheme for second order derivative at grid point `I`
along dimension `dim`. E.g. if `dim=1`, this approximates `∂ₓₓ`.
"""
function D2⁻⁻(ϕ,I,dim,h=step(ϕ,dim))
    Im   = decrement_index(I,dim,1)
    Imm  = decrement_index(I,dim,2)
    return (ϕ[Imm] - 2ϕ[Im] + ϕ[I]) / h^2
end

"""
    curvature(ϕ,I)

Compute ∇ ⋅ (∇ ϕ / |∇ ϕ|) at grid point `I`.
"""
function curvature(ϕ,I)
    N = ambient_dimension(ϕ)
    if N == 2
        ϕx  = D⁰(ϕ,I,1)
        ϕy  = D⁰(ϕ,I,2)
        ϕxx = D2⁰(ϕ,I,1)
        ϕyy = D2⁰(ϕ,I,2)
        ϕxy = D2(ϕ,I,(2,1))
        κ   = (ϕxx*(ϕy)^2 - 2*ϕy*ϕx*ϕxy + ϕyy*ϕx^2) / (ϕx^2 + ϕy^2)^(3/2)
        return κ
    elseif N == 3
        ϕx  = D⁰(ϕ,I,1)
        ϕy  = D⁰(ϕ,I,2)
        ϕz  = D⁰(ϕ,I,3)
        ϕxx = D2⁰(ϕ,I,1)
        ϕyy = D2⁰(ϕ,I,2)
        ϕzz = D2⁰(ϕ,I,3)
        ϕxy = D2(ϕ,I,(2,1))
        ϕxz = D2(ϕ,I,(3,1))
        # TODO: test + simplify this
        κ   = (ϕxx*(ϕy)^2 - 2*ϕy*ϕx*ϕxy + ϕyy*ϕx^2 + ϕx^2*ϕzz - 2*ϕx*ϕz*ϕxz + ϕz^2*ϕxx + ϕy^2*ϕzz - 2*ϕy*ϕz*ϕyz + ϕz^2*ϕyy) / (ϕx^2 + ϕy^2)^3/2
        return κ
    else
        notimplemented()
    end
end
