"""
    struct LevelSetEquation

Representation of a level-set equation of the form `ϕₜ + ∑ᵢ (Fᵢ) = 0`, where
each `Fᵢ` is a [`LevelSetTerm`](@ref).

A `LevelSetEquation` has a `current_state` representing a level-set function at
the `current_time`. It can be stepped foward in time using
`integrate!(ls,Δt)`, which evolves the level set equation for a time interval `Δt`,
modifying in the process its `current_state` and `current_time`.

Boundary conditions are specified in the field `bc`, and the scheme for the
time-integration can be set in the `integrator` field. Finally, the `cfl` fields
can be used to control the multiplier in the `cfl` condition; that is, `Δt = cfl*Δt_max`.
"""
Base.@kwdef mutable struct LevelSetEquation
    terms::Tuple
    integrator::TimeIntegrator
    levelset::LevelSet
    t::Float64 = 0
    buffers = allocate_buffers(levelset,integrator)
    cfl::Float64 = 0.5
    boundary_condition::BoundaryCondition
end

function allocate_buffers(ϕ::LevelSet,::ForwardEuler)
    f = levelset_function(ϕ)
    return (deepcopy(f),)
end

function allocate_buffers(ϕ::LevelSet,::RK2)
    f = levelset_function(ϕ)
    return (deepcopy(f),deepcopy(f))
end

function Base.show(io::IO, eq::LevelSetEquation)
    print(io, "Level-set equation given by\n")
    print(io, "\n \t ϕₜ + ")
    terms = eq.terms
    for term in terms[1:end-1]
        print(io, term)
        print(io, " + ")
    end
    print(io, terms[end])
    print(io, " = 0")
    print(io, "\n\n Current time $(eq.t)")
    return io
end

# getters
levelset(eq::LevelSetEquation) = eq.levelset
levelset_function(eq::LevelSetEquation) = levelset_function(eq.levelset)
vals(eq::LevelSetEquation) = eq |> levelset_function |> vals
mesh(eq::LevelSetEquation) = eq |> levelset_function |> mesh
current_time(eq::LevelSetEquation) = eq.t
buffers(eq::LevelSetEquation) = eq.buffers
time_integrator(eq::LevelSetEquation) = eq.integrator
terms(eq::LevelSetEquation) = eq.terms
cfl(eq::LevelSetEquation) = eq.cfl
boundary_condition(eq::LevelSetEquation) = eq.boundary_condition

function compute_rhs!(eq::LevelSetEquation)
    ϕ  = levelset_function(eq)
    bc = boundary_condition(eq)
    dϕ = buffers(eq)[1]
    t  = current_time(eq)
    pars = (bc=boundary_condition(eq), terms=terms(eq))
    compute_rhs!(dϕ, ϕ, pars, t) # write in dϕ
    return eq
end

function compute_rhs!(dϕ, ϕ, p, t)
    applybc!(ϕ, p.bc)
    nodes = ϕ |> vals_mesh |> NodeIterator
    for I in interior_indices(ϕ, p.bc)
        x = nodes[I]
        _compute_terms(p.terms, ϕ, I)
        dϕ[I] = -_compute_terms(p.terms, ϕ, I)
    end
    return dϕ
end

function _compute_terms(terms::Tuple, ϕ, I)
    sum(terms) do term
        _compute_term(term, ϕ, I)
    end
end

"""
    integrate!(eq::LevelSetEquation,tf,Δt=Inf)

Integrate the [`LevelSetEquation`](@ref) `ls` up to time `tf`,
mutating the `current_state` and `current_time` of the object `eq` in the
process.

An optional parameter `Δt` can be passed to specify a maximum time-step
allowed for the integration. Note that the internal time-steps taken to evolve
the level-set up to `tf` may be smaller than `Δt` due to stability reasons
related to the `terms` and `integrator` field in `ls`.
"""
function integrate!(eq::LevelSetEquation, tf, Δt=Inf)
    tc = current_time(eq)
    msg = "final time $(tf) must be larger than the initial time $(tc):
           the level-set equation cannot be solved back in time"
    @assert tf >= tc msg
    ϕ = levelset_function(eq)
    pars = (bc=boundary_condition(eq), terms=terms(eq), mesh=mesh(eq),
        cfl=cfl(eq), tc=tc, tf=tf, Δt=Δt, buffers=buffers(eq))
    # dynamic dispatch. Should not be a problem provided enough computation is
    # done inside the function below
    ϕ, tc = _integrate!(ϕ, pars, time_integrator(eq))
    eq.t = tc
    #
    return eq
end

function _integrate!(ϕ, p, ::ForwardEuler)
    Δt_cfl = p.cfl * compute_cfl(p.terms, ϕ, p.bc)
    Δt = min(p.Δt, Δt_cfl) # minimum of desired and required time step
    tc = p.tc
    dϕ = p.buffers[1]
    while tc <= p.tf - eps(tc)
        Δt = min(Δt, p.tf - tc) # if needed, take a smaller time-step to exactly land on tf
        applybc!(ϕ, p.bc)
        compute_rhs!(dϕ, ϕ, p, tc)
        for I in interior_indices(ϕ, p.bc)
            ϕ[I] = ϕ[I] + Δt * dϕ[I]
        end
        tc += Δt
        @debug tc, Δt
    end
    return ϕ, tc
end

function _integrate!(ϕ, p, ::RK2)
    Δt_cfl = p.cfl * compute_cfl(p.terms, ϕ, p.bc)
    Δt = min(p.Δt, Δt_cfl)
    tc = p.tc
    dϕ = p.buffers[1]
    ϕmid = p.buffers[2]
    while tc <= p.tf - eps(tc)
        Δt = min(Δt, p.tf - tc) # if needed, take a smaller time-step to exactly land on tf
        applybc!(ϕ, p.bc)
        compute_rhs!(dϕ, ϕ, p, tc)
        for I in interior_indices(ϕ, p.bc)
            ϕmid[I] = ϕ[I] + 0.5 * Δt * dϕ[I]
        end
        applybc!(ϕmid, p.bc)
        compute_rhs!(dϕ, ϕmid, p, tc+0.5*Δt)
        for I in interior_indices(ϕ, p.bc)
            ϕ[I] = ϕ[I] + Δt * dϕ[I]
        end
        tc += Δt
        @debug tc, Δt
    end
    return ϕ, tc
end

"""
    abstract type LevelSetTerm

A typical term in a level-set evolution equation.
"""
abstract type LevelSetTerm end

function compute_cfl(terms, ϕ, bc::BoundaryCondition)
    minimum(terms) do term
        _compute_cfl(term, ϕ, bc::BoundaryCondition)
    end
end

# generic method, loops over indices
function _compute_cfl(term::LevelSetTerm, ϕ, bc::BoundaryCondition)
    dt = Inf
    iter = ϕ |> mesh |> NodeIterator
    for I in interior_indices(iter, bc)
        cfl = _compute_cfl(term, ϕ, I)
        dt = min(dt, cfl)
    end
    return dt
end

# generic method, loops over dimensions
function _compute_cfl(term::LevelSetTerm, ϕ, I)
    N = ambient_dimension(ϕ)
    minimum(1:N) do dim
        _compute_cfl(term, ϕ, I, dim)
    end
end

"""
    struct AdvectionTerm{V,M} <: LevelSetTerm

Level-set advection term representing  `𝐯 ⋅ ∇ϕ`.
"""
Base.@kwdef struct AdvectionTerm{F,S<:SpatialScheme} <: LevelSetTerm
    velocity::F
    scheme::S = Upwind()
end
velocity(adv::AdvectionTerm) = adv.velocity
scheme(adv::AdvectionTerm) = adv.scheme

Base.show(io::IO, ::AdvectionTerm) = print(io, "u⃗ ⋅ ∇ ϕ")

function _compute_term(term::AdvectionTerm, ϕ::CartesianGridFunction, I)
    N = ambient_dimension(ϕ)
    u = velocity(term)[I]
    sum(1:N) do dim
        _compute_term(term, ϕ, I, dim, u[dim])
    end
end

@inline function _compute_term(term::AdvectionTerm, ϕ::CartesianGridFunction, I, dim, v)
    sch = scheme(term)
    # for dimension dim, compute the upwind derivative and multiply by the
    # velocity v in that dimension
    if sch === Upwind()
        if v > 0
            return v * D⁻(ϕ, I, dim)
        else
            return v * D⁺(ϕ, I, dim)
        end
    elseif sch === WENO5()
        if v > 0
            return v * weno5⁻(ϕ, I, dim)
        else
            return v * weno5⁺(ϕ, I, dim)
        end
    else
        error("scheme $sch not implemented")
    end
end


function _compute_cfl(term::AdvectionTerm, ϕ, I, dim)
    𝐮 = velocity(term)[I]
    # for each dimension, compute the upwind derivative and multiply by the
    # velocity and add to buffer
    Δx = step(ϕ)[dim]
    return Δx / abs(𝐮[dim])
end

"""
    struct CurvatureTerm{V,M} <: LevelSetTerm

Level-set curvature term representing `bκ|∇ϕ|`, where `κ = ∇ ⋅ (∇ϕ/|∇ϕ|) ` is
the curvature.
"""
struct CurvatureTerm{N,T,V} <: LevelSetTerm
    b::CartesianGridFunction{N,T,V}
end
coefficient(cterm::CurvatureTerm) = cterm.b

Base.show(io::IO, ::CurvatureTerm) = print(io, "b κ|∇ϕ|")

function _compute_term(term::CurvatureTerm, ϕ, I)
    N = ambient_dimension(ϕ)
    b = coefficient(term)
    κ = curvature(ϕ, I)
    # compute |∇ϕ|
    ϕ2 = sum(1:N) do dim
        D⁰(ϕ, I, dim)^2
    end
    # update
    return b[I] * κ * sqrt(ϕ2)
end

function _compute_cfl(term::CurvatureTerm, ϕ, I, dim)
    b = coefficient(term)[I]
    Δx = step(ϕ)[dim]
    return (Δx)^2 / (2 * abs(b))
end

"""
    struct NormalMotionTerm{V,M} <: LevelSetTerm

Level-set advection term representing  `v |∇ϕ|`. This `LevelSetTerm` should be
used for internally generated velocity fields; for externally generated
velocities you may use [`AdvectionTerm`](@ref) instead.
"""
Base.@kwdef struct NormalMotionTerm{T} <: LevelSetTerm
    speed::T
end
speed(adv::NormalMotionTerm) = adv.speed

Base.show(io::IO, ::NormalMotionTerm) = print(io, "v|∇ϕ|")

function _compute_term(term::NormalMotionTerm, ϕ::CartesianGridFunction, I)
    v = speed(term)[I]
    ∇ = _compute_∇_norm(v, ϕ, I)
    return ∇ * v
end

function _compute_∇_norm(v, ϕ, I)
    N = ambient_dimension(ϕ)
    mA0², mB0² = sum(1:N) do dim
        h = step(ϕ, dim)
        A = D⁻(ϕ, I, dim) + 0.5 * h * limiter(D2⁻⁻(ϕ, I, dim), D2⁰(ϕ, I, dim))
        B = D⁺(ϕ, I, dim) - 0.5 * h * limiter(D2⁺⁺(ϕ, I, dim), D2⁰(ϕ, I, dim))
        if v > 0.0
            SVector(positive(A)^2, negative(B)^2)
        else
            SVector(negative(A)^2, positive(B)^2)
        end
    end
    return sqrt(mA0² + mB0²)
end

function _compute_cfl(term::NormalMotionTerm, ϕ, I, dim)
    u = speed(term)[I]
    Δx = step(ϕ)[dim]
    return Δx / abs(u)
end

@inline positive(x) = max(x, zero(x))
@inline negative(x) = min(x, zero(x))

# eq. (6.28)
function limiter(x, y)
    x * y < zero(x) || return zero(x)
    return abs(x) <= abs(y) ? x : y
end

"""
    struct ReinitializationTerm <: LevelSetTerm

Level-set term representing  `sign(ϕ) (|∇ϕ| - 1)`. This `LevelSetTerm` should be
used for reinitializing the level set into a signed distance function: for a
sufficiently large number of time steps this term allows one to solve the
Eikonal equation |∇ϕ| = 1.
"""
Base.@kwdef struct ReinitializationTerm <: LevelSetTerm
end

Base.show(io::IO, ::ReinitializationTerm) = print(io, "sign(ϕ) (|∇ϕ| - 1)")

function _compute_term(term::ReinitializationTerm, ϕ::CartesianGridFunction, I)
    v = sign(ϕ[I])
    ∇ = _compute_∇_norm(v, ϕ, I)
    return (∇ - 1.0) * v
end

_compute_cfl(::ReinitializationTerm, ϕ, I, dim) = step(ϕ)[dim]
