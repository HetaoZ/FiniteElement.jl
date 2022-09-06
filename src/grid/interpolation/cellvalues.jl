function CellScalarValues(quad_rule::QuadratureRule, func_interpol::Interpolation,
    geom_interpol::Interpolation=func_interpol)
CellScalarValues(Float64, quad_rule, func_interpol, geom_interpol)
end

function CellScalarValues(::Type{T}, quad_rule::QuadratureRule{dim,shape}, func_interpol::Interpolation,
    geom_interpol::Interpolation=func_interpol) where {dim,T,shape<:AbstractRefShape}

@assert getdim(func_interpol) == getdim(geom_interpol)
@assert getrefshape(func_interpol) == getrefshape(geom_interpol) == shape
n_qpoints = length(getweights(quad_rule))

# Function interpolation
n_func_basefuncs = getnbasefunctions(func_interpol)
N    = fill(zero(T)          * T(NaN), n_func_basefuncs, n_qpoints)
dNdx = fill(zero(Vec{dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)
dNdξ = fill(zero(Vec{dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)

# Geometry interpolation
n_geom_basefuncs = getnbasefunctions(geom_interpol)
M    = fill(zero(T)          * T(NaN), n_geom_basefuncs, n_qpoints)
dMdξ = fill(zero(Vec{dim,T}) * T(NaN), n_geom_basefuncs, n_qpoints)

for (qp, ξ) in enumerate(quad_rule.points)
    for i in 1:n_func_basefuncs
        dNdξ[i, qp], N[i, qp] = gradient(ξ -> value(func_interpol, i, ξ), ξ, :all)
    end
    for i in 1:n_geom_basefuncs
        dMdξ[i, qp], M[i, qp] = gradient(ξ -> value(geom_interpol, i, ξ), ξ, :all)
    end
end

detJdV = fill(T(NaN), n_qpoints)

CellScalarValues{dim,T,shape}(N, dNdx, dNdξ, detJdV, M, dMdξ, quad_rule, func_interpol, geom_interpol)
end

# CellVectorValues
struct CellVectorValues{dim,T<:Real,refshape<:AbstractRefShape,M} <: CellValues{dim,T,refshape}
N::Matrix{Vec{dim,T}}
dNdx::Matrix{Tensor{2,dim,T,M}}
dNdξ::Matrix{Tensor{2,dim,T,M}}
detJdV::Vector{T}
M::Matrix{T}
dMdξ::Matrix{Vec{dim,T}}
qr::QuadratureRule{dim,refshape,T}
# The following fields are deliberately abstract -- they are never used in
# performance critical code, just stored here for convenience.
func_interp::Interpolation{dim,refshape}
geo_interp::Interpolation{dim,refshape}
end

function CellVectorValues(quad_rule::QuadratureRule, func_interpol::Interpolation, geom_interpol::Interpolation=func_interpol)
CellVectorValues(Float64, quad_rule, func_interpol, geom_interpol)
end

function CellVectorValues(::Type{T}, quad_rule::QuadratureRule{dim,shape}, func_interpol::Interpolation,
    geom_interpol::Interpolation=func_interpol) where {dim,T,shape<:AbstractRefShape}

@assert getdim(func_interpol) == getdim(geom_interpol)
@assert getrefshape(func_interpol) == getrefshape(geom_interpol) == shape
n_qpoints = length(getweights(quad_rule))

# Function interpolation
n_func_basefuncs = getnbasefunctions(func_interpol) * dim
N    = fill(zero(Vec{dim,T})      * T(NaN), n_func_basefuncs, n_qpoints)
dNdx = fill(zero(Tensor{2,dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)
dNdξ = fill(zero(Tensor{2,dim,T}) * T(NaN), n_func_basefuncs, n_qpoints)

# Geometry interpolation
n_geom_basefuncs = getnbasefunctions(geom_interpol)
M    = fill(zero(T)          * T(NaN), n_geom_basefuncs, n_qpoints)
dMdξ = fill(zero(Vec{dim,T}) * T(NaN), n_geom_basefuncs, n_qpoints)

for (qp, ξ) in enumerate(quad_rule.points)
    basefunc_count = 1
    for basefunc in 1:getnbasefunctions(func_interpol)
        dNdξ_temp, N_temp = gradient(ξ -> value(func_interpol, basefunc, ξ), ξ, :all)
        for comp in 1:dim
            N_comp = zeros(T, dim)
            N_comp[comp] = N_temp
            N[basefunc_count, qp] = Vec{dim,T}((N_comp...,))

            dN_comp = zeros(T, dim, dim)
            dN_comp[comp, :] = dNdξ_temp
            dNdξ[basefunc_count, qp] = Tensor{2,dim,T}((dN_comp...,))
            basefunc_count += 1
        end
    end
    for basefunc in 1:n_geom_basefuncs
        dMdξ[basefunc, qp], M[basefunc, qp] = gradient(ξ -> value(geom_interpol, basefunc, ξ), ξ, :all)
    end
end

detJdV = fill(T(NaN), n_qpoints)
MM = Tensors.n_components(Tensors.get_base(eltype(dNdx)))

CellVectorValues{dim,T,shape,MM}(N, dNdx, dNdξ, detJdV, M, dMdξ, quad_rule, func_interpol, geom_interpol)
end

function reinit!(cv::CellValues{dim}, x::AbstractVector{Vec{dim,T}}) where {dim,T}
n_geom_basefuncs = getngeobasefunctions(cv)
n_func_basefuncs = getnbasefunctions(cv)
length(x) == n_geom_basefuncs || throw_incompatible_coord_length(length(x), n_geom_basefuncs)

@inbounds for i in 1:length(cv.qr.weights)
    w = cv.qr.weights[i]
    fecv_J = zero(Tensor{2,dim})
    for j in 1:n_geom_basefuncs
        fecv_J += x[j] ⊗ cv.dMdξ[j, i]
    end
    detJ = det(fecv_J)
    detJ > 0.0 || throw_detJ_not_pos(detJ)
    cv.detJdV[i] = detJ * w
    Jinv = inv(fecv_J)
    for j in 1:n_func_basefuncs
        cv.dNdx[j, i] = cv.dNdξ[j, i] ⋅ Jinv
    end
end
end

getngeobasefunctions(cv::Values) = size(cv.M, 1)

"""
    getnquadpoints(fe_v::Values)

Return the number of quadrature points for the `Values` object.
"""
getnquadpoints(fe::Values) = length(fe.qr.weights)

"""
    reinit!(cv::CellValues, x::Vector)
    reinit!(bv::FaceValues, x::Vector, face::Int)

Update the `CellValues`/`FaceValues` object for a cell or face with coordinates `x`.
The derivatives of the shape functions, and the new integration weights are computed.
"""
function reinit!(cv::CellValues{dim}, x::AbstractVector{Vec{dim,T}}) where {dim,T}
    n_geom_basefuncs = getngeobasefunctions(cv)
    n_func_basefuncs = getnbasefunctions(cv)
    length(x) == n_geom_basefuncs || throw_incompatible_coord_length(length(x), n_geom_basefuncs)

    @inbounds for i in 1:length(cv.qr.weights)
        w = cv.qr.weights[i]
        fecv_J = zero(Tensor{2,dim})
        for j in 1:n_geom_basefuncs
            fecv_J += x[j] ⊗ cv.dMdξ[j, i]
        end
        detJ = det(fecv_J)
        detJ > 0.0 || throw_detJ_not_pos(detJ)
        cv.detJdV[i] = detJ * w
        Jinv = inv(fecv_J)
        for j in 1:n_func_basefuncs
            cv.dNdx[j, i] = cv.dNdξ[j, i] ⋅ Jinv
        end
    end
end