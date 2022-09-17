# 预处理只需要创建网格
" 1/2/3维  Line2/Quad4/Hex8 单元网格"
function generate(g::RectangularGrid{dim,T}) where dim where T <: CubeElement
    nnode = g.nel .+ 1
    d = (g.stop .- g.start) ./ g.nel 
    coords = ntuple(axis -> [g.start[axis] + d[axis] * (i-1) for i = 1:nnode[axis]], dim)

    nodes = Array{Node{dim}}(undef, nnode)
    for c in CartesianIndices(nodes)
        node_id = cartesian_to_node_id(c, nnode)
        x0 = Vec(ntuple(i -> coords[i][c[i]], dim))
        x = Vec(ntuple(i -> coords[i][c[i]], dim))
        nodes[c] = Node{dim}(node_id, x0, x, tensorzeros(dim),tensorzeros(dim),tensorzeros(dim),tensorzeros(dim))
    end

    elements = Array{T}(undef,g.nel)
    for c in CartesianIndices(elements)
        if T == Line
            node1 = nodes[c].id
            node2 = nodes[c[1]+1].id
            connection = (node1, node2)
            faces = ((node1,), (node2,))
        elseif T == Quadrilateral
            node1 = nodes[c].id
            node2 = nodes[c[1]+1, c[2]].id
            node3 = nodes[c[1]+1, c[2]+1].id
            node4 = nodes[c[1], c[2]+1].id
            connection = (node1, node2, node3, node4)
            faces = ((node1,node2), (node2,node3), (node3,node4), (node4,node1))
        elseif T == Hexahedron
            P1 = nodes[c].id
            P2 = nodes[c[1]+1, c[2], c[3]].id
            P3 = nodes[c[1]+1, c[2]+1, c[3]].id
            P4 = nodes[c[1], c[2]+1, c[3]].id

            P5 = nodes[c[1], c[2], c[3]+1].id
            P6 = nodes[c[1]+1, c[2], c[3]+1].id
            P7 = nodes[c[1]+1, c[2]+1, c[3]+1].id
            P8 = nodes[c[1], c[2]+1, c[3]+1].id  
            
            connection = (P1,P2,P3,P4,P5,P6,P7,P8)
            faces = ((P1,P5,P8,P4), 
                     (P2,P3,P7,P6), 
                     (P1,P2,P6,P5), 
                     (P3,P4,P8,P7), 
                     (P1,P4,P3,P2), 
                     (P5,P6,P7,P8))
        else
            error("undefined element type")
        end
        

        elements[c] = Element(T, connection, faces)
    end

    nodes = reshape(nodes, (length(nodes),))
    elements = reshape(elements, (length(elements),))

    grid = Grid{dim,T}(nodes, elements, SurfaceTopology(T), g)
    get_surface_topo!(grid)

    return grid
end

@inline cartesian_to_node_id(c::CartesianIndex{3}, nnode) = (c[3]-1) * nnode[1]*nnode[2] + (c[2]-1)*nnode[1] + c[1]
@inline cartesian_to_node_id(c::CartesianIndex{2}, nnode) = (c[2]-1)*nnode[1] + c[1]
@inline cartesian_to_node_id(c::CartesianIndex{1}, nnode) = c[1]

"获取指定空间范围内包含的结点编号向量"
function find_nodes(s::Structure, start, stop)
    return find_nodes(s.grid, start, stop)
end

function find_nodes(s::Structure, P1, P2, P3)
    return find_nodes(s.grid, P1, P2, P3)
end


# ---------------------------------------------------
# 创建Structure实例

"创建 Structure 实例，可根据 grid_prototype 自动生成网格"
function Structure(material::AbstractMaterial, grid_prototype::AbstractGridPrototype, solver::AbstractSolver)
    return Structure(material, generate(grid_prototype), solver)
end

"创建 Structure 实例"
function Structure(material::AbstractMaterial, grid::Grid{dim,T}, solver::AbstractSolver) where {dim,T}
    
    s = Structure{dim}(material, grid, solver, new_states(material, getnelems(grid), getnq(grid)), new_solution(solver, getndofs(grid)), AbstractConstrain[], Dict(), true)

    update_states!(s)
    return s
end