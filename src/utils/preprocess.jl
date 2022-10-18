# 预处理只需要创建网格
" 1/2/3维  Line2/Quad4/Hex8 单元网格"
function generate(g::RectangularGrid{dim,T}) where dim where T <: CubeElement
    
    nodes = generate_nodes(g)
    elements = generate_cube_elements(g, nodes)

    nodes = reshape(nodes, (length(nodes),))
    elements = reshape(elements, (length(elements),))

    grid = Grid{dim,T}(nodes, elements, SurfaceTopology(T), g)
    get_surface_topo!(grid)

    return grid
end

" 2/3维  Tri3/Tet4 单元网格"
function generate(g::RectangularGrid{dim,T}) where dim where T <: SimplexElement
    
    nodes = generate_nodes(g)
    elements = generate_simplex_elements(g, nodes)

    nodes = reshape(nodes, (length(nodes),))
    elements = reshape(elements, (length(elements),))

    grid = Grid{dim,T}(nodes, elements, SurfaceTopology(T), g)
    get_surface_topo!(grid)

    return grid
end

function generate_nodes(g::RectangularGrid{dim,T}) where dim where T
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
    return nodes
end

function generate_cube_elements(g::RectangularGrid{dim,T}, nodes) where dim where T 
    elements = Array{TO_CUBE_TYPE[T]}(undef,g.nel)
    for c in CartesianIndices(elements)
        if T == Line
            node1 = nodes[c].id
            node2 = nodes[c[1]+1].id
            connection = (node1, node2)
            faces = ((node1,), (node2,))
        elseif T in (Quadrilateral, Triangle)
            node1 = nodes[c].id
            node2 = nodes[c[1]+1, c[2]].id
            node3 = nodes[c[1]+1, c[2]+1].id
            node4 = nodes[c[1], c[2]+1].id
            connection = (node1, node2, node3, node4)
            faces = ((node1,node2), (node2,node3), (node3,node4), (node4,node1))
        elseif T in (Hexahedron, Tetrahedron)
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
        
        elements[c] = Element(TO_CUBE_TYPE[T], connection, faces)
    end
    return elements
end

const TO_CUBE_TYPE = Dict(Triangle=>Quadrilateral, Quadrilateral=>Quadrilateral, Tetrahedron=>Hexahedron, Hexahedron=>Hexahedron)
const N_SIMPLEX = (1,2,5) # number of simplexes divided from a k-dimensional cube
const TRI = ([1,2,3], [2,4,3]) # Tri's connections om a rect
const TRI_FACE = (([1,2],[2,3],[3,1]), ([2,4],[4,3],[3,1]))
const TET = ([1,2,3,5], [2,4,3,8], [6,8,2,5], [7,5,3,8], [5,8,2,3]) # Tet's connections in a cube
const TET_FACE = ntuple(i->([TET[i][1],TET[i][3],TET[i][2]], 
[TET[i][1],TET[i][2],TET[i][4]],
[TET[i][1],TET[i][4],TET[i][3]],
[TET[i][2],TET[i][3],TET[i][4]]), length(TET))

function generate_simplex_elements(g::RectangularGrid{dim,T}, nodes) where dim where T <: SimplexElement

    elements = generate_cube_elements(g, nodes)
    simplexes = Array{T}(undef,N_SIMPLEX[dim],size(elements)...)

    for c in CartesianIndices(elements)
        if T == Triangle
            P =  [   nodes[c             ].id,
                     nodes[c[1]+1, c[2]  ].id,
                     nodes[c[1],   c[2]+1].id,
                     nodes[c[1]+1, c[2]+1].id]

            for j in 1:2
                simplexes[j,c] = Element(T, Tuple(P[TRI[j]]), ntuple(i->Tuple(P[TRI_FACE[j][i]]), 3))
            end
        elseif T == Tetrahedron
            P = [
                nodes[c                     ].id,
                nodes[c[1]+1, c[2],   c[3]  ].id,
                nodes[c[1]  , c[2]+1, c[3]  ].id,
                nodes[c[1]+1, c[2]+1, c[3]  ].id,
                nodes[c[1],   c[2],   c[3]+1].id,
                nodes[c[1]+1, c[2],   c[3]+1].id,
                nodes[c[1]  , c[2]+1, c[3]+1].id,
                nodes[c[1]+1, c[2]+1, c[3]+1].id
                ]
            
            for j in 1:5
                simplexes[j,c] = Element(T, Tuple(P[TET[j]]), ntuple(i->Tuple(P[TET_FACE[j][i]]), 4))
            end
        else
            println("elemtype is ", T)
            error("undefined element type")
        end
    end

    return simplexes
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