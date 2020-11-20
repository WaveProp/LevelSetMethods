abstract type LevelSetTerm{N} end

dimension(::LevelSetTerm{N}) where {N} = N


struct Advection{N} <: LevelSetTerm{N}
    v
    method
end  

dimension(adv::Advection{N}) where {N} = N+1

V   = #champs de vitesse 
adv = Advection(V,:upwind)

adv(ϕ)

struct LevelSet{A,B}
    ϕ::A
    grid::B
end

function upwind(ϕ::LevelSet{A,CartesianGrid})

end    

function upwind(ϕ::Level{A,Octree})
        
end 

abstract type Mesh end

function (adv::Advection)(buffer,ϕ::LevelSet)
    if adv.method == :upwind
        upwind!(buffer,ϕ,adv.v)        
    end        
    return buffer
end    

struct NormalAvection
    a #scalar
    met
end    

function compute_rhs(ϕ,terms::Tuple)
    buffer = zero(ϕ)
    for term in terms
        term(buffer,ϕ)
    end 
    return buffer       
end    









