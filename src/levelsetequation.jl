"""
    struct LevelSetEquation

Representation of a level-set equation of the form `ϕₜ + ∑ᵢ (Fᵢ) = 0`, where
each `Fᵢ` is a `LevelSetTerm`.

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
    state::LevelSet
    t::Float64=0
    # TODO: buffer allocation depends on the `integrator` field, with some
    # integrators requiring e.g. two copies of `state` for the  intermediate
    # stages. Maybe write a `preallocate_buffer` function for handling that
    # automatically in the constructor?
    buffer=zero(state)
    cfl::Float64=0.5
end

function Base.show(io::IO, eq::LevelSetEquation)
    print(io,"Level-set equation given by\n")
    print(io,"\n \t ϕₜ + ")
    terms = eq.terms
    for term in terms[1:end-1]
        print(io,term)
        print(io," + ")
    end
    print(io,terms[end])
    print(io," = 0")
    print(io,"\n\n Current time $(eq.t)")
    return io
end

# getters
current_state(ls::LevelSetEquation) = ls.state
current_time(ls::LevelSetEquation)  = ls.t[]
buffer(ls::LevelSetEquation)        = ls.buffer
time_integrator(ls::LevelSetEquation) = ls.integrator
terms(ls::LevelSetEquation)           = ls.terms
cfl(ls::LevelSetEquation)  = ls.cfl

"""
    integrate!(ls::LevelSetEquation,tf,Δt=Inf)

Integrate the [`LevelSetEquation`](@ref) `ls` up to time `tf`,
mutating the `current_state` and `current_time` of the object `ls` in the
process.

An optional parameter `Δt` can be passed to specify a maximum time-step
allowed for the integration. Note that the internal time-steps taken to evolve
the level-set up to `tf` may be smaller than `Δt` due to stability reasons
related to the `terms` and `integrator` field in `ls`.
"""
function integrate!(ls::LevelSetEquation,tf,Δt=Inf)
    tc = current_time(ls)
    α  = cfl(ls)
    msg = "final time $(tf) must be larger than the initial time $(tc):
           the level-set equation cannot be solved back in time"
    @assert tf >= tc msg
    b          = buffer(ls)
    ϕ          = current_state(ls)
    integrator = time_integrator(ls)

    # dynamic dispatch. Should not be a problem provided enough computation is
    # done inside the function below
    ϕ,b = _integrate!(ϕ,b,integrator,α,ls.terms,tc,tf,Δt)
    ls.t = tf
    # reassigning ls.ϕ and ls.b needed because these could have been
    # swapped inside of _integrate!
    ls.state  = ϕ
    ls.buffer = b
    return ls
end

function _integrate!(ϕ::LevelSet,buffer::LevelSet,integrator::ForwardEuler,α,terms,tc,tf,Δt)
    Δt_cfl = α * compute_cfl(terms,ϕ)
    Δt     = min(Δt,Δt_cfl)
    while tc <= tf - eps(tc)
        Δt  = min(Δt,tf-tc) # if needed, take a smaller time-step to exactly land on tf
        applybc!(ϕ)
        for I in interior_indices(ϕ)
            buffer[I] = _compute_terms(terms,ϕ,I)
            buffer[I] = ϕ[I] - Δt * buffer[I] # muladd?
        end
        ϕ, buffer = buffer,ϕ # swap the roles, no copies
        tc += Δt
        @debug tc,Δt
    end
    # @assert tc ≈ tf
    return ϕ,buffer
end

function _integrate!(ϕ::LevelSet,buffers,integrator::RK2,α,terms,tc,tf,Δt)
    buffer1,buffer2 = buffers[1],buffers[2]
    Δt_cfl = α * compute_cfl(terms,ϕ)
    Δt     = min(Δt,Δt_cfl)
    while tc <= tf - eps(tc)
        Δt  = min(Δt,tf-tc) # if needed, take a smaller time-step to exactly land on tf
        applybc!(ϕ)
        m    = mesh(ϕ)
        iter = NodeIterator(m)
        for I in interior_indices(ϕ)
            tmp = _compute_terms(terms,ϕ,I)
            buffer1[I] = ϕ[I] - Δt * tmp # muladd?
            buffer2[I] = ϕ[I] - 0.5*Δt * tmp # muladd?
        end
        applybc!(buffer1)
        for I in interior_indices(ϕ)
            tmp = _compute_terms(terms,buffer1,I)
            buffer2[I] -= 0.5*Δt * tmp
        end
        ϕ,buffer1,buffer2 = buffer2,ϕ,buffer1 # swap the roles, no copies
        tc += Δt
        @debug tc,Δt
    end
    # @assert tc ≈ tf
    return ϕ,(buffer1,buffer2)
end

function _integrate!(ϕ::LevelSet,buffer::LevelSet,integrator::RKLM2,α,terms,tc,tf,Δt)
    Δt_cfl = α * compute_cfl(terms,ϕ)
    Δt     = min(Δt,Δt_cfl)
    while tc <= tf - eps(tc)
        Δt  = min(Δt,tf-tc) # if needed, take a smaller time-step to exactly land on tf
        applybc!(ϕ)
        grid = mesh(ϕ)
        for I in interior_indices(ϕ)
            tmp = _compute_terms(terms,ϕ,I)
            buffer[I] = tmp
        end
        for I in interior_indices(ϕ)
            ϕ[I] = ϕ[I] - Δt * buffer[I] # muladd?
        end
        applybc!(ϕ)
        for I in interior_indices(ϕ)
            tmp = _compute_terms(terms,ϕ,I)
            buffer[I] = ϕ[I] - 0.5*Δt*tmp + 0.5*Δt*buffer[I]
        end
        ϕ,buffer = buffer,ϕ # swap the roles, no copies
        tc += Δt
        @debug tc,Δt
    end
    # @assert tc ≈ tf
    return ϕ,buffer
end

# for a tuple of terms, sum their contributions
function _compute_terms(terms::Tuple,ϕ,I)
    sum(terms) do term
        _compute_term(term,ϕ,I)
    end
end

# recipes for Plots
@recipe function f(eq::LevelSetEquation)
    ϕ = current_state(eq)
    t = current_time(eq)
    N = ambient_dimension(ϕ)
    if N == 2 # 2d contour plot
        seriestype --> :contour
        levels --> [0,]
        aspect_ratio --> :equal
        colorbar --> false
        title --> "t = $t"
        # seriescolor --> :black
        return ϕ
    else
        notimplemented()
    end
end
