function upwind(ϕ,u,grid::CartesianGrid)
    setbc!(ϕ,:periodic1)    
    h = gridsize(grid)    
    T = eltype(ϕ)    
    m,n = size(grid)
    ∇ϕ  = zeros(SVector{2,Float64},m,n)
    for i in 2:m-1
        for j in 2:n-1
            # ϕ\_x
            if u[i,j][1]>0
                ϕx = (ϕ[i,j] - ϕ[i-1,j]) / h[2]
            else
                ϕx = (ϕ[i+1,j] - ϕ[i,j]) / h[2]
            end        
            # ϕy
            if u[i,j][2]>0
                ϕy  = (ϕ[i,j] - ϕ[i,j-1]) / h[2]
            else
                ϕy = (ϕ[i,j+1] - ϕ[i,j]) / h[2]
            end        
            ∇ϕ[i,j] = (ϕx,ϕy)
        end
    end
    return ∇ϕ
end    

function setbc!(ϕ,bctype)
    if bctype == :periodic1
        ϕ[1,:]   = ϕ[end-1,:]
        ϕ[end,:] = ϕ[2,:]
        ϕ[:,1]   = ϕ[:,end-1]
        ϕ[:,end] = ϕ[:,2]
    elseif bctype == :neumann1
        ϕ[1,:]   = ϕ[2,:]
        ϕ[end,:] = ϕ[end-1,:]
        ϕ[:,1]   = ϕ[:,2]
        ϕ[:,end] = ϕ[:,end-1]
    else 
        error("uknown boundary condition $bctype")    
    end    
    return ϕ
end    

