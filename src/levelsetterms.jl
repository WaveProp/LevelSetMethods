"""
    abstract type LevelSetTerm

A typical term in a level-set evolution equation.
"""
abstract type LevelSetTerm end

"""
    compute_terms(terms,ϕ,bc)
    compute_terms!(buffer,terms,ϕ,bc)

Given a tuple `terms` containing `LevelSetTerm`s, compute the contribution of all
these terms to the level set equation. A `buffer` can be passed to allocate the output.
"""
function compute_terms!(buffer::MeshField,terms::Tuple,ϕ::LevelSet)
    applybc!(ϕ) # TODO: who is actually responsible for calling this?
    grid = mesh(ϕ)
    for I in interior_indices(ϕ)
        buffer[I] = sum(terms) do term
            _compute_term!(term,ϕ,I)    
        end            
    end   
    return buffer
end    
compute_terms(args...) = compute_terms!(zero(ϕ),args...)

"""
    struct AdvectionTerm{V,M} <: LevelSetTerm

Level-set advection term representing  `𝐯 ⋅ ∇ϕ`.
"""
Base.@kwdef struct AdvectionTerm{V,M} <: LevelSetTerm
    velocity::MeshField{V,M}
end
velocity(adv::AdvectionTerm) = adv.velocity

function _compute_term(term::AdvectionTerm,ϕ,I,dim)
    𝐮 = velocity(term)
    N = dimension(ϕ)
    # for dimension dim, compute the upwind derivative and multiply by the
    # velocity
    v = 𝐮[I][dim]
    if v > 0
        return v*D⁻(ϕ,I,dim)
    else
        return v*D⁺(ϕ,I,dim)
    end
end

function _compute_term(term::AdvectionTerm,ϕ,I)
    N = dimension(term)    
    sum(1:N) do dim
        _compute_term(term,ϕ,I,dim)    
    end    
end

function _compute_cfl(term::AdvectionTerm,ϕ,I,dim)
    𝐮 = velocity(term)[I]
    N = dimension(ϕ)
    # for each dimension, compute the upwind derivative and multiply by the
    # velocity and add to buffer
    Δx = meshsize(ϕ)[dim]
    return Δx/abs(𝐮[dim])
end    

function _compute_cfl(term::AdvectionTerm,ϕ,I)
    N = dimension(term)    
    minimum(1:N) do dim
        _compute_cfl(term,ϕ,I,dim)    
    end    
end    

# generic method, loops over dimensions
function _compute_cfl(term::LevelSetTerm,ϕ::LevelSet,I)
    N = dimension(ϕ)    
    minimum(1:N) do dim
        _compute_cfl(term,ϕ,I,dim)
    end
end

# generic method, loops over indices
function _compute_cfl(term::LevelSetTerm,ϕ::LevelSet)
    minimum(interior_indices(ϕ)) do I
        _compute_cfl(term,ϕ,I)    
    end    
end    

"""
    struct CurvatureTerm{V,M} <: LevelSetTerm

Level-set curvature term representing `bκ|∇ϕ|`, where `κ = ∇ ⋅ (∇ϕ/|∇ϕ|) ` is
the curvature.
"""
struct CurvatureTerm{V,M} <: LevelSetTerm
    b::MeshField{V,M}
end
coefficient(cterm::CurvatureTerm) = cterm.b

function _compute_term(term::CurvatureTerm,ϕ,I)
    b = coefficient(term)
    κ = curvature(ϕ,I)
    # compute |∇ϕ|
    ∇ϕ = map(1:N) do dim
        D⁰(ϕ,I,dim)
    end
    # update
    buffer[I] += b[I]*κ*norm(∇ϕ,2)
    return buffer
end

function _compute_cfl(term::CurvatureTerm,ϕ,I,dim)
    b = coefficient(term)[I]
    Δx = meshsize(ϕ)[dim]
    return (Δx)^2/(2*abs(b))
end    

function curvature(ϕ::LevelSet,I)
    N = dimension(ϕ)
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

"""
    struct NormalAdvectionTerm{V,M} <: LevelSetTerm

Level-set advection term representing  `v |∇ϕ|`. This `LevelSetTerm` should be
used for internally generated velocity fields; for externally generated
velocities you may use `AdvectionTerm` instead.
"""
@Base.kwdef struct NormalAdvectionTerm{V,M} <: LevelSetTerm
    speed::MeshField{V,M}
end
speed(adv::NormalAdvectionTerm) = adv.speed

function _compute_term(term::NormalAdvectionTerm,ϕ,I)
    u = speed(term)
    N = dimension(ϕ)
    v = u[I]
    mA0² = 0.0
    mB0² = 0.0
    for dim in 1:N
        h = meshsize(ϕ,dim)
        # eq. (6.22-6.27) generalized for any dimensions
        A = D⁻(ϕ,I,dim) + 0.5 * h * limiter(D2⁻⁻(ϕ,I,dim), D2⁰(ϕ,I,dim))
        B = D⁺(ϕ,I,dim) - 0.5 * h * limiter(D2⁺⁺(ϕ,I,dim), D2⁰(ϕ,I,dim))
        if v > 0.0
            mA0² += positive(A)^2
            mB0² += negative(B)^2
        else
            mA0² += negative(A)^2
            mB0² += positive(B)^2
        end
    end
    ∇ = sqrt(mA0² + mB0²)
    return ∇ * v
end

function _compute_cfl(term::NormalAdvectionTerm,ϕ,I,dim)
    u = speed(term)[I]
    Δx = meshsize(ϕ)[dim]
    return Δx/abs(u) 
end

@inline positive(x) = x > zero(x) ? x : zero(x)
@inline negative(x) = x < zero(x) ? x : zero(x)

# eq. (6.20-6.21)
function g(x, y)
    tmp = zero(x)
    if x > zero(x); tmp += x*x; end
    if y < zero(x); tmp += y*y; end
    return sqrt(tmp)
end

# eq. (6.28)
function limiter(x, y)
    x*y < zero(x) || return zero(x)
    return abs(x) <= abs(y) ? x : y
end
