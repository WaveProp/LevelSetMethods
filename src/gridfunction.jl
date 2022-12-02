"""
    struct CartesianGridFunction{N,T,V}

A function `f : ℝᴺ → V` defined by discrete values on the nodes of a
`UniformCartesianMesh`. The `order` parameter specifies the polynomial order of
the interpolation scheme.
"""
struct CartesianGridFunction{N,T,V}
    vals::Array{V,N}
    vals_mesh::UniformCartesianMesh{N,T} # the mesh where `vals` are given
    els_mesh::UniformCartesianMesh{N,T}  # mesh for the elements2
    order::Int
end

mesh(f::CartesianGridFunction) = f.els_mesh
vals(f::CartesianGridFunction) = f.vals
meshsize(f::CartesianGridFunction, args...) = step(mesh(f), args...)
ambient_dimension(f::CartesianGridFunction{N}) where {N} = N

domain(f::CartesianGridFunction) = domain(mesh(f))

function CartesianGridFunction(f::Function, U::HyperRectangle{N,T}; meshsize,
                               order=1) where {N,T}
    els_mesh = UniformCartesianMesh(U; step=meshsize)
    sz = size(els_mesh)
    vals_mesh = UniformCartesianMesh(U, sz .* order)
    # resize mesh so that it is divisible by order
    vals = [f(x) for x in nodes(vals_mesh)]
    return CartesianGridFunction(vals, vals_mesh, els_mesh, order)
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
    return (bernstein_interpolant(f, I, vandermond) for I in CartesianIndices(els))
end
