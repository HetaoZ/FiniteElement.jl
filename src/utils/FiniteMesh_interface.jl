using FiniteMesh 

"*.msh file format supported"
function read_grid(msh_file::String)
    mesh = Mesh(msh_file)
    grid = mesh_to_grid(mesh)
    return grid
end

function mesh_to_grid(mesh)
    dim = size(mesh.points,2)
    nodes = [Node(dim,i,Vec(Tuple(mesh.points[i,:])...)) for i in axes(mesh.points,1)]
    nel, nen = size(mesh.cells.index[end])
    elemtype = infer_elemtype(dim, nen)
    elements = [Element(elemtype, Tuple(mesh.cells.index[end][i,:])) for i in 1:nel]
    return Grid{dim,elemtype}(nodes, elements)
end

const DIM_ELEMTYPE = Dict(
    (1,2) => Line,
    (2,3) => Triangle,
    (2,4) => Quadrilateral,
    (3,4) => Tetrahedron,
    (3,8) => Hexahedron
)

function infer_elemtype(dim, nen)
    return DIM_ELEMTYPE[(dim, nen)]
end