# 单元类型
include("element/element.jl")
# 母单元上的高斯积分
include("quadrature/quadrature.jl")
# 母单元上的插值,和母单元到物理单元的映射?
include("interpolation/interpolation.jl")

Node(dim::Int) = Node{dim}(0,tensorzeros(dim),tensorzeros(dim),tensorzeros(dim),tensorzeros(dim),tensorzeros(dim),tensorzeros(dim))

getdim(node::Node{dim}) where dim = dim
getdofs(node::Node{dim}) where dim = collect((node.id-1)*getdim(node)+1:node.id*getdim(node))

getdim(::Grid{dim,T}) where {dim,T} = dim
getnnodes(g::Grid) = length(g.nodes)
getnelems(g::Grid) = length(g.elements)
getndofs(g::Grid{dim,T}) where {dim,T} = dim * length(g.nodes)


"获取指定空间范围内包含的结点编号向量"
function find_nodes(g::Grid, start, stop)
    nodeids = Int[]
    for node in g.nodes
        if between(node.x, start, stop)
            push!(nodeids, node.id)
        end
    end
    return nodeids
end

