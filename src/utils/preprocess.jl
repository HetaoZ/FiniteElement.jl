# 预处理只需要创建网格

function generate(g::RectangularGrid{2,Quadrilateral})
    dim = 2
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

    elements = Array{Quadrilateral}(undef,g.nel)
    for c in CartesianIndices(elements)
        node1 = nodes[c].id
        node2 = nodes[c[1]+1, c[2]].id
        node3 = nodes[c[1]+1, c[2]+1].id
        node4 = nodes[c[1], c[2]+1].id
        connection = (node1, node2, node3, node4)
        faces = ((node1,node2), (node2,node3), (node3,node4), (node4,node1))
        elements[c] = Element(Quadrilateral, connection, faces)
    end

    nodes = reshape(nodes, (length(nodes),))
    elements = reshape(elements, (length(elements),))
    return Grid{2,Quadrilateral}(nodes, elements)
end

function cartesian_to_node_id(c::CartesianIndex{3}, nnode)
    return (c[3]-1) * nnode[1]*nnode[2] + (c[2]-1)*nnode[1] + c[1]
end

function cartesian_to_node_id(c::CartesianIndex{2}, nnode)
    return (c[2]-1)*nnode[1] + c[1]
end

function cartesian_to_node_id(c::CartesianIndex{1}, nnode)
    return c[1]
end

function Structure(material::AbstractMaterial, grid::Grid{dim,T}, solver::AbstractSolver) where {dim,T}
    return Structure{dim}(material, grid, solver, new_states(material, getnelems(grid), getnq(grid)), new_solution(solver, getndofs(grid)), AbstractConstrain[], Dict(), true)
end

"根据prototype自动生成网格"
function Structure(material::AbstractMaterial, grid_prototype::AbstractGridPrototype, solver::AbstractSolver)
    return Structure(material, generate(grid_prototype), solver)
end

"获取指定空间范围内包含的结点编号向量"
function find_nodes(s::Structure, start, stop)
    return find_nodes(s.grid, start, stop)
end

