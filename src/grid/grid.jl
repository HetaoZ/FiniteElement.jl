# 单元类型
include("element/element.jl")
# 母单元上的高斯积分
include("quadrature/quadrature.jl")
# 母单元上的插值,和母单元到物理单元的映射?
include("interpolation/interpolation.jl")

"不使用原型网格， 已知 surface"
function Grid{dim,T}(nodes, elements, surface_topo) where dim where T <: AbstractElementType
    return Grid{dim,T}(nodes, elements, surface_topo, NonePrototype{dim}())
end

"不使用原型网格，且暂无定义 surface "
function Grid{dim,T}(nodes, elements) where dim where T <: AbstractElementType
    surface_topo = SurfaceTopology(eltype(elements))
    return Grid{dim,T}(nodes, elements, surface_topo, NonePrototype{dim}())
end

SurfaceTopology(::Type{Element{dim,N,M,L}}) where {dim,N,M,L} = SurfaceTopology(L)
SurfaceTopology(L::Int) = SurfaceTopology(zeros(Int,L,0))

Node(dim::Int) = Node{dim}(0,tensorzeros(dim),tensorzeros(dim),tensorzeros(dim),tensorzeros(dim),tensorzeros(dim),tensorzeros(dim))

Node(dim::Int,id,x0) = Node{dim}(id,x0,x0,tensorzeros(dim),tensorzeros(dim),tensorzeros(dim),tensorzeros(dim))

getdim(node::Node{dim}) where dim = dim
getdofs(node::Node{dim}) where dim = collect((node.id-1)*getdim(node)+1:node.id*getdim(node))

getdim(::Grid{dim,T}) where {dim,T} = dim
getnnodes(g::Grid) = length(g.nodes)
getnelems(g::Grid) = length(g.elements)
getndofs(g::Grid{dim,T}) where {dim,T} = dim * length(g.nodes)
getnq(g::Grid) = length(g.elements[1].quad_rule.weights)


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

"对于三维空间，可根据三点确定的平面来搜寻平面附近的结点"
function find_nodes(g::Grid{3,T}, P1, P2, P3) where T<:AbstractElementType
    nodeids = Int[]
    for node in g.nodes
        if point_near_plane(node.x, P1,P2,P3)
            push!(nodeids, node.id)
        end
    end
    return nodeids
end

"允许对 grid 进行平移、旋转、缩放"
function transform!(grid::Grid{dim,T}; translation = (0,0,0), rotation = ((0,0,1), 0), scaling = ((0,0,0), (1,1,1))) where {dim,T}
    translate_grid!(grid, translation)
    rotate_grid!(grid, rotation[1], rotation[2])
    zoom_grid!(grid, scaling[1], scaling[2])
end

function vector2vec(v)
    return Vec(Tuple(v)...)
end

function translate_grid!(grid::Grid{dim,T}, v) where {dim,T}
    for i in eachindex(grid.nodes)
        grid.nodes[i].x0 = vector2vec(grid.nodes[i].x0 .+ v)
        grid.nodes[i].x = vector2vec(grid.nodes[i].x .+ v)
    end
    if typeof(grid.prototype) <: RectangularGrid
        start = start .+ v
        stop = stop .+ v
        grid.prototype = RectangularGrid{dim,T}(start, stop, grid.prototype.nel)
    end
end

function quaternion(t)
    return Quaternion(0, expand(t))
end

function expand(t)
    v = zeros(Float64,3)
    for i in eachindex(t)
        v[i] = t[i]
    end
    return v
end

function rotate_grid!(grid::Grid{dim,T}, rotation_axis::NTuple{3,R}, rotation_angle::Real) where {dim,T,R<:Real}
    # 四元数旋转得到角度是2φ，所以除以2
    φ = rotation_angle * 0.5

    # 旋转四元数
    q = Quaternion(cos(φ), sin(φ) * collect(rotation_axis))
    q_conj = conj(q)
    # 旋转函数
    q_rotate(x) = imag_part(q * quaternion(x) * q_conj)[1:length(x)]

    for i in eachindex(grid.nodes)
        grid.nodes[i].x0 = vector2vec(q_rotate(grid.nodes[i].x0))
        grid.nodes[i].x = vector2vec(q_rotate(grid.nodes[i].x))
    end
    if typeof(grid.prototype) <: RectangularGrid
        start = q_rotate(start)
        stop = q_rotate(stop)
        grid.prototype = RectangularGrid{dim,T}(start, stop, grid.prototype.nel)
    end
end

"以 origin 为基点的缩放操作"
function zoom_grid!(grid::Grid{dim,T}, origin, scale) where {dim,T}
    for i in eachindex(grid.nodes)
        grid.nodes[i].x0 = vector2vec(zoom(grid.nodes[i].x0, origin, scale))
        grid.nodes[i].x = vector2vec(zoom(grid.nodes[i].x, origin, scale))
        grid.nodes[i].d = vector2vec(zoom(grid.nodes[i].d, origin, scale))
        grid.nodes[i].u = vector2vec(zoom(grid.nodes[i].u, origin, scale))
        grid.nodes[i].a = vector2vec(zoom(grid.nodes[i].a, origin, scale))
    end
    if typeof(grid.prototype) <: RectangularGrid
        start = zoom(grid.prototype.start, origin, scale)
        stop = zoom(grid.prototype.stop, origin, scale)
        grid.prototype = RectangularGrid{dim,T}(start, stop, grid.prototype.nel)
    end
end

function zoom(x, origin, scale)
    return (x .- origin) .* scale .+ origin
end