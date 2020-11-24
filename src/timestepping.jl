abstract type TimeIntegrator end

struct ForwardEuler
end

function evolve!(buffer,ϕ,terms,bc::BoundaryCondition,dt)
    compute_terms!(buffer,terms,ϕ,bc)
    @. ϕ.vals = ϕ.vals - dt*buffer.vals
    return ϕ
end


