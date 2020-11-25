abstract type TimeIntegrator end

Base.@kwdef struct ForwardEuler{T}
    cfl::Float64 = 0.5
    buffer::T
end

function evolve!(ϕ,feuler::ForwardEuler,terms,bc,t,Δtmax=Inf)
    buffer      = feuler.buffer
    fill!(values(buffer),0)
    buffer, Δtˢ = compute_terms!(buffer,terms,ϕ,bc)    
    α = feuler.cfl
    Δt = min(Δtmax,α*Δtˢ)
    # buffer now holds the rhs. Update ϕ in-place
    axpy!(-Δt,buffer.vals,ϕ.vals) # ϕ = ϕ - dt*buffer
    return ϕ,t+Δt
end    

Base.@kwdef struct RK2{T}
    cfl::Float64 = 0.5 
    buffers::NTuple{2,T}
end

function evolve!(ϕ,integ::RK2,terms,bc,t,Δtmax=Inf)
    α = integ.cfl    
    buffer1,buffer2 = integ.buffers[1],integ.buffers[2]
    fill!(values(buffer1),0)
    fill!(values(buffer2),0)
    #
    buffer1, Δtˢ = compute_terms!(buffer1,terms,ϕ,bc) 
    Δt = min(Δtmax,α*Δtˢ) 
    axpy!(-Δt/2,buffer1.vals,ϕ.vals) # ϕ = ϕ - dt/2*buffer1
    #
    @. buffer1.vals = ϕ.vals + Δt * buffer1.vals   
    buffer2, _   = compute_terms!(buffer2,terms,buffer1,bc)    
    #    
    axpy!(-Δt/2,buffer2.vals,ϕ.vals) # ϕ = ϕ - dt/2*buffer2
    return ϕ,t+Δt
end    


