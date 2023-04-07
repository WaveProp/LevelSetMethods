"""
    struct CartesianGridFunction{N,T,V}

A piecewise-polynomial function `p : ℝᴺ → V` defined by discrete values on the
nodes of a `UniformCartesianMesh`. The `order` parameter specifies the
polynomial order of the interpolation scheme.
"""
struct CartesianGridFunction{N,T,V}
    vals::Array{V,N}
    vals_mesh::UniformCartesianMesh{N,T} # the mesh where `vals` are given
    els_mesh::UniformCartesianMesh{N,T}  # mesh for the elements2
    interpolants::Dict{CartesianIndex{N}, BernsteinPolynomial{N,V}}
    order::Int
end

mesh(f::CartesianGridFunction) = f.els_mesh
vals(f::CartesianGridFunction) = f.vals
vals_mesh(f::CartesianGridFunction) = f.vals_mesh
interpolants(f::CartesianGridFunction) = f.interpolants
interpolant(f::CartesianGridFunction, I::CartesianIndex) = f.interpolants[I] 
ambient_dimension(f::CartesianGridFunction{N}) where {N} = N

Base.size(f::CartesianGridFunction) = size(vals(f))
Base.getindex(f::CartesianGridFunction, args...) = getindex(vals(f), args...)
Base.setindex!(f::CartesianGridFunction, v, args...) = setindex!(vals(f), v, args...)

Base.step(f::CartesianGridFunction) = step(vals_mesh(f))
Base.step(f::CartesianGridFunction,dim) = step(vals_mesh(f),dim)

domain(f::CartesianGridFunction) = domain(mesh(f))

"""
    CartesianGridFunction(f::Function, U::HyperRectangle{N,T}; meshsize,
    order=1)

Project the function `f` on a `UniformCartesianGrid` of `U` with meshsize given
by `meshsize`. The `order` parameter specifies the polynomial order of the
interpolation scheme for reconstructing a continuous function from the grid
values.
"""
function CartesianGridFunction(f::Function, U::HyperRectangle{N,T}; meshsize,
                               order=1) where {N,T}
    els_mesh = UniformCartesianMesh(U; step=meshsize)
    sz = size(els_mesh)
    vals_mesh = UniformCartesianMesh(U, sz .* order)
    # resize mesh so that it is divisible by order
    vals = [f(x) for x in nodes(vals_mesh)]

    # vandermond matrix
    nodes1d = collect(range(0, 1, order + 1))
    buffer = zeros(Float64, ntuple(i -> order + 1, N)) #FIXME typage
    vandermond = ones(float(T), length(buffer), length(buffer))
    idxs = CartesianIndices(buffer)
    for (i, I) in enumerate(idxs)
        for (j, J) in enumerate(idxs)
            for d in 1:N
                n = order
                k = J[d] - 1
                xd = nodes1d[I[d]]
                vandermond[i, j] *= binomial(n, k) * xd^k * (1 - xd)^(n - k)
            end
        end
    end
    # factor the matrix
    F = lu!(vandermond)

    els = ElementIterator(els_mesh)
    interpolants = Dict(I => bernstein_interpolant(vals, order, I, els[I], F) for I in CartesianIndices(els))

    return CartesianGridFunction(vals, vals_mesh, els_mesh, interpolants, order)
end

function bernstein_interpolant(vals::Array{V,N}, order::Int, I::CartesianIndex{N}, U, vandermond) where {V,N}
    coeffs = zeros(V, ntuple(i -> order + 1, N))
    idxs_vals = CartesianIndices(ntuple(N) do d
            return ((I[d] - 1) * order + 1):(I[d] * order + 1)
        end)
    @views ldiv!(vec(coeffs), vandermond, vec(vals[idxs_vals]))
    return BernsteinPolynomial(coeffs, U)
end

function merge_interpolant(f::CartesianGridFunction{N,T,V}, I::CartesianIndex{N}, J::CartesianIndex{N}) where {N,T,V}
    ax = findfirst(x->!iszero(x), Tuple(I) .- Tuple(J))
    # vandermond matrix
    nodes1d = collect(range(0, 1, f.order + 1)); nodes1d_ax = collect(range(0, 1, f.order*2 + 1))
    buffer = zeros(V, ntuple(i -> i == ax ? f.order*2 + 1 : f.order + 1, N))
    vandermond = ones(float(T), length(buffer), length(buffer))
    idxs = CartesianIndices(buffer)
    for (i, I) in enumerate(idxs)
        for (j, J) in enumerate(idxs)
            for d in 1:N
                n  = d == ax ? f.order*2 : f.order
                k  = J[d] - 1
                xd = d == ax ? nodes1d_ax[I[d]] : nodes1d[I[d]]
                vandermond[i, j] *= binomial(n, k) * xd^k * (1 - xd)^(n - k)
            end
        end
    end
    # factor the matrix
    F = lu!(vandermond)
    
    coeffs = zeros(V, ntuple(i -> i == ax ? f.order*2 + 1 : f.order + 1, N))
    idxs_vals = CartesianIndices(ntuple(N) do d
            if d == ax
                l = min(I[d], J[d]); u = max(I[d], J[d])
                return ((l - 1) * f.order + 1):(u * f.order + 1)                    
            end
            return ((I[d] - 1) * f.order + 1):(I[d] * f.order + 1)
        end)
    @views ldiv!(vec(coeffs), F, vec(vals(f)[idxs_vals]))

    els = ElementIterator(f.els_mesh)
    U1 = els[I]; U2 = els[J]
    U  = HyperRectangle(min.(low_corner(U1), low_corner(U2)), max.(high_corner(U1), high_corner(U2)))
    return BernsteinPolynomial(coeffs, U)
end

function bernstein_interpolant(f::CartesianGridFunction{N,T,V}, I::CartesianIndex{N},
                               vandermond=bernstein_interpolation_matrix(f)) where {N,T,V}
    coeffs = zeros(V, ntuple(i -> f.order + 1, N))
    els = ElementIterator(f.els_mesh)
    U = els[I]
    idxs_vals = CartesianIndices(ntuple(N) do d
                                     return ((I[d] - 1) * f.order + 1):(I[d] * f.order + 1)
                                 end)
    @views ldiv!(vec(coeffs), vandermond, vec(vals(f)[idxs_vals]))
    return BernsteinPolynomial(coeffs, U)
end

function bernstein_interpolation_matrix(f::CartesianGridFunction{N,T,V}) where {N,T,V}
    nodes1d = collect(range(0, 1, f.order + 1))
    buffer = zeros(V, ntuple(i -> f.order + 1, N))
    # vandermond matrix
    vandermond = ones(float(T), length(buffer), length(buffer))
    idxs = CartesianIndices(buffer)
    for (i, I) in enumerate(idxs)
        for (j, J) in enumerate(idxs)
            for d in 1:N
                n = f.order
                k = J[d] - 1
                xd = nodes1d[I[d]]
                vandermond[i, j] *= binomial(n, k) * xd^k * (1 - xd)^(n - k)
            end
        end
    end
    # factor the matrix
    F = lu!(vandermond)
    return F
end

function bernstein_interpolants(f::CartesianGridFunction{N,T,V}) where {N,T,V}
    vandermond = bernstein_interpolation_matrix(f)
    els = ElementIterator(f.els_mesh)
    return [bernstein_interpolant(f, I, vandermond) for I in CartesianIndices(els)]
end

element_index_for_point(x::SVector{N}, f::CartesianGridFunction{N}) where {N} = element_index_for_point(x, mesh(f))

function curvature(f::CartesianGridFunction{N,T,V}, x::SVector{N}) where {N,T,V}
    I = element_index_for_point(x, f)
    p = bernstein_interpolant(f, I)
    return curvature(p, x)
end
