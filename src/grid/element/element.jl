# 具体定义各类单元的函数
include("quadrilateral.jl")


# 一般单元类型的函数
function fetch_data(elem::Element{dim,N,M,L},nodes::Vector{Node},field) where {dim,N,M,L}
    data = zeros(Float64,dim,N)
    for j in 1:N
        data[:,j] = getfield(nodes[elem.connection[j]], field)
    end
    return data
end

elem_x0(elem::Element, nodes::Vector{Node}) = fetch_data(elem, nodes, :x0)
elem_x(elem::Element, nodes::Vector{Node}) = fetch_data(elem, nodes, :x)
elem_d(elem::Element, nodes::Vector{Node}) = fetch_data(elem, nodes, :d)

function getdofs(elem::Element{dim,N,M,L}) where {dim,N,M,L}
    dofs = zeros(Int,dim,N)
    for i in 1:dim, j in 1:N
        dofs[i,j] = (elem.connection[j]-1) * dim + i
    end
    return collect(dofs)
end

"""
膨胀系数等于当前体积除以初始体积
"""
function expansion_coeff(elem::Element, nodes)
    v0 = init_volume(elem, nodes)
    v = volume(elem, nodes)
    return v/v0
end

