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
Base.step(f::CartesianGridFunction, args...) = step(mesh(f), args...)
ambient_dimension(f::CartesianGridFunction{N}) where {N} = N

domain(f::CartesianGridFunction) = domain(mesh(f))

function CartesianGridFunction(f::Function, U::HyperRectangle{N,T}; step,
                               order=1) where {N,T}
    vals_mesh = UniformCartesianMesh(U; step)
    # resize mesh so that it is divisible by order
    sz = ceil.(Int, size(vals_mesh) ./ order) .* order
    vals_mesh = UniformCartesianMesh(U, sz)
    vals = [f(x) for x in nodes(vals_mesh)]
    els_mesh = UniformCartesianMesh(U, div.(sz, order))
    return CartesianGridFunction(vals, vals_mesh, els_mesh, order)
end

# monomial interpolants
function monomial_interpolant!(p̂,Fin,f::CartesianGridFunction{N,T,V}, I::CartesianIndex{N}) where {N,T,V}
    p̂,F = preallocate_monomial_interpolant(f)
    @assert F == Fin
    c = coeffs(p̂)
    els = ElementIterator(f.els_mesh)
    U   = els[I]
    idxs = CartesianIndices(ntuple(N) do d
                                return ((I[d] - 1) * f.order + 1):(I[d] * f.order + 1)
                            end)
    ldiv!(vec(c),F,vec(vals(f)[idxs]))
    p̂  = power2bernstein(c)
    p  = BernsteinPolynomial(p̂.coeffs,p̂.degree,U)
    return U,p
end

function monomial_interpolant(f::CartesianGridFunction,I::CartesianIndex)
    p̂,F = preallocate_monomial_interpolant(f)
    monomial_interpolant!(p̂,F,f,I)
end

function preallocate_monomial_interpolant(f::CartesianGridFunction{N,T,V}) where {N,T,V}
    nodes1d = collect(range(0, 1, f.order + 1))
    p̂       = Polynomial(zeros(V, ntuple(i -> f.order + 1, N)))
    # vandermond matrix
    c = coeffs(p̂)
    vandermond = ones(float(T), length(c), length(c))
    idxs = CartesianIndices(c)
    for (i,I) in enumerate(idxs)
        for (j,J) in enumerate(idxs)
            for d in 1:N
                vandermond[i,j] *= nodes1d[I[d]]^(J[d]-1)
            end
        end
    end
    # factor the matrix
    F = lu(vandermond)
    return p̂,F
end

function monomial_interpolants(f::CartesianGridFunction{N,T,V}) where {N,T,V}
    p̂,F = preallocate_monomial_interpolant(f)
    els = ElementIterator(f.els_mesh)
    return (monomial_interpolant!(p̂, F, f, I) for I in CartesianIndices(els))
end


# Lagrange interpolants
function preallocate_lagrange_interpolant(f::CartesianGridFunction{N,T,V}) where {N,T,V}
    nodes1d = ntuple(i -> collect(range(0, 1, f.order + 1)), N)
    p̂ = TensorLagInterp(zeros(V, ntuple(i -> f.order + 1, N)), nodes1d)
    return p̂
end

function lagrange_interpolant!(p̂::TensorLagInterp{N}, f::CartesianGridFunction{N},
                      I::CartesianIndex{N}) where {N}
    p = f.order
    els = ElementIterator(f.els_mesh)
    U = els[I]
    idxs = CartesianIndices(ntuple(N) do d
                                return ((I[d] - 1) * p + 1):(I[d] * p + 1)
                            end)
    copy!(p̂.vals, f.vals[idxs])
    xl = low_corner(U)
    w = width(U)
    return U, (x) -> p̂((x - xl) ./ w)
end
function lagrange_interpolant(f::CartesianGridFunction{N}, I::CartesianIndex{N}) where {N}
    p̂ = preallocate_lagrange_interpolant(f)
    return lagrange_interpolant!(p̂, f, I)
end

function lagrange_interpolants(f::CartesianGridFunction{N,T,V}) where {N,T,V}
    nodes1d = ntuple(i -> collect(range(0, 1, f.order + 1)), N)
    p̂ = TensorLagInterp(zeros(V, ntuple(i -> f.order + 1, N)), nodes1d)
    els = ElementIterator(f.els_mesh)
    return (lagrange_interpolant!(p̂, f, I) for I in CartesianIndices(els))
end
