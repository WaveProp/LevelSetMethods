"""
    struct LevelSetEquation

Representation of a level-set equation of the form `ϕₜ + sum(terms) = 0`, where 
each `t ∈ terms` is a `LevelSetTerm`. 

A `LevelSetEquation` has a `current_state` representing a level-set function at
the `current_time`. It can also be stepped foward in time using
`evolve!(ls,Δt)`, which evolves the level set equation for a time interval `Δt`,
modifying in the process its `current_state` and `current_time`. 

Boundary conditions are specified in the field `bc`, and the scheme for the
time-integration can be set in the `integrator` field.
"""
Base.@kwdef mutable struct LevelSetEquation
    terms::Tuple
    integrator::TimeIntegrator
    state::LevelSet
    t::Float64
    buffer
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
time_integrator(ls)                 = ls.integrator
terms(ls)                           = ls.terms

"""
    integrate!(ls::LevelSetEquation,tf,Δt=Inf)

Integrate the [`LevelSetEquation`](@ref) `ls` up to time `tf`, 
mutating the `current_state` and `current_time` of the object `ls` in the
process. 

An optional parameter `Δt` can be passed to specify a maximum time-step
allowed for the integration. Note that the internal time-steps taken to evolve
the level-set up to `tf` may be smaller than `Δt` due to stability reasons
related to the `terms` and `integrator` employed.
"""
function integrate!(ls::LevelSetEquation,tf,Δt=Inf)
    tc = current_time(ls)
    msg = "final time $(tf) must be larger than the initial time $(tc): 
           the level-set equation cannot be solved back in time"
    @assert tf >= tc msg
    b          = buffer(ls)
    ϕ          = current_state(ls)
    integrator = time_integrator(ls)

    # dynamic dispatch. Should not be a problem provided enough computation is
    # done inside of the function below
    ϕ,b = _integrate!(ϕ,b,integrator,ls.terms,tc,tf,Δt)
    ls.t = tf
    # reassigning ls.ϕ and ls.b may be needed because these could have been
    # swapped inside of _integrate!
    ls.state  = ϕ
    ls.buffer = b
    return ls
end

function _integrate!(ϕ::LevelSet,buffer::LevelSet,integrator::ForwardEuler,terms,tc,tf,Δt)
    α      = cfl(integrator)
    Δt_cfl = α * compute_cfl(terms,ϕ) 
    Δt     = min(Δt,Δt_cfl)
    while tc <= tf - eps(tc)
        Δt  = min(Δt,tf-tc) # if needed, take a smaller time-step to exactly land on tf        
        applybc!(ϕ) 
        grid = mesh(ϕ)
        for I in interior_indices(ϕ)
            buffer[I] = sum(terms) do term
                _compute_term(term,ϕ,I)    
            end
            buffer[I] = ϕ[I] - Δt * buffer[I] # muladd?
        end       
        ϕ, buffer = buffer,ϕ # swap the roles, no copies
        tc += Δt
        @debug tc,Δt
    end
    # @assert tc ≈ tf
    return ϕ,buffer
end

# recipes for Plots
@recipe function f(eq::LevelSetEquation)
    ϕ = current_state(eq)    
    t = current_time(eq)
    N = dimension(ϕ)
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