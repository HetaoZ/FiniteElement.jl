# 具体定义各类单元的函数
include("line2.jl")
include("tri3.jl")
include("quad4.jl")
include("tet4.jl")
include("hex8.jl")


# 一般单元类型的函数
"""
返回 (dim,N) 形状的数据
"""
function fetch_data(elem::Element{dim,N,M,L},nodes::Vector{Node{dim}},field) where {dim,N,M,L}
    data = zeros(Float64,dim,N)
    for j in 1:N
        data[:,j] = getfield(nodes[elem.connection[j]], field)
    end
    return data
end

elem_x0(elem::Element{dim}, nodes::Vector{Node{dim}}) where dim = fetch_data(elem, nodes, :x0)
elem_x(elem::Element{dim}, nodes::Vector{Node{dim}}) where dim = fetch_data(elem, nodes, :x)
elem_d(elem::Element{dim}, nodes::Vector{Node{dim}}) where dim = fetch_data(elem, nodes, :d)

# function getcoords(elem::Element{dim,N,M,L}, nodes::Vector{Node{dim}}) where {dim,N,M,L}
#     x = Vector{Vec{dim,Float64}}(undef,length(elem.connection))
#     for (i, node_id) in enumerate(elem.connection)
#         x[i] = nodes[node_id].x
#     end
#     return x
# end

function getdofs(elem::Element{dim,N,M,L}) where {dim,N,M,L}
    dofs = zeros(Int,dim,N)
    for i in 1:dim, j in 1:N
        dofs[i,j] = (elem.connection[j]-1) * dim + i
    end
    return vec(dofs)
end

"""
膨胀系数等于当前体积除以初始体积
"""
function expansion_coeff(elem::Element, nodes)
    v0 = init_volume(elem, nodes)
    v = volume(elem, nodes)
    return v/v0
end

function elem_density(elem::Element, nodes, material::AbstractMaterial)
    return material.ρ₀ / expansion_coeff(elem, nodes)
end

# -----------------------------
# CellValues 
function CellScalarValues(quad_rule::QuadratureRule{dim,shape}, ip::Interpolation) where {dim,shape<:AbstractRefShape}

    # Function interpolation
    n = getnbasefunctions(ip) # 插值点数（=结点数）
    nq = length(getweights(quad_rule))

    N    = fill(zero(Float64) * Float64(NaN), n, nq)  # size = (n,nq)
    dNdx = fill(zero(Float64) * Float64(NaN), dim, n, nq)  # size = (dim,n,nq)
    dNdξ = fill(zero(Float64) * Float64(NaN), dim, n, nq)  # size = (dim,n,nq)

    for (qp, ξ) in enumerate(quad_rule.points)
        for i in 1:n
            dNdξ[:, i, qp], N[i, qp] = gradient(ξ -> value(ip, i, ξ), ξ, :all)
        end
    end

    detJdV = fill(Float64(NaN), nq)

    CellScalarValues{dim,Float64,shape}(N, dNdx, dNdξ, detJdV)
end

"""
    reinit!(cv::CellValues, ip::Interpolation, quad_rule::QuadratureRule, x::Matrix)

Update the `CellValues` object for a cell with coordinates `x`, where size(x) = (n, dim).
The derivatives of the shape functions, and the new integration weights are computed.

[@ref] 王勖成《有限单元法》P133
"""
function reinit!(cv::CellValues{dim}, quad_rule::QuadratureRule, x::Matrix{Float64}) where {dim}

    nq = length(quad_rule.weights) # 积分点数

    @inbounds for i_qpoint in 1:nq
        # 计算第i个积分点上的 detJ * dV
        # size(J) = (dim,dim)
        J = cv.dNdξ[:,:,i_qpoint] * x
        detJ = det(J)
        detJ > 0.0 || throw_detJ_not_pos(detJ)
        cv.detJdV[i_qpoint] = detJ * quad_rule.weights[i_qpoint]

        Jinv = inv(J)
        cv.dNdx[:,:,i_qpoint] = Jinv * cv.dNdξ[:,:,i_qpoint]

        # println()
        # println("qpoint = ", quad_rule.points[i_qpoint])
        # println("dNdξ = ")
        # display(cv.dNdξ[:,:,i_qpoint]); println()
        # println("detJ = ", det(J))
        
    end
end

@noinline throw_detJ_not_pos(detJ) = throw(ArgumentError("det(J) is not positive: det(J) = $(detJ)"))

"""
d 代表单元各结点位移，size(d) = (n,dim)。与之对应，size(cv.dNdx) = (dim,n,nq)
"""
function disp_gradient(cv::CellValues{dim}, q_point::Int, d::Matrix{Float64}) where {dim}
    ∇d = Tensor{2,dim,Float64}(shape_gradient(cv, q_point) * d)
    return ∇d
end

@propagate_inbounds shape_gradient(cv::CellValues, q_point::Int) = cv.dNdx[:,:, q_point]

@propagate_inbounds getdetJdV(cv::CellValues, q_point::Int) = cv.detJdV[q_point]

@propagate_inbounds shape_value(cv::CellValues, i::Int, q_point::Int) = cv.N[i, q_point]