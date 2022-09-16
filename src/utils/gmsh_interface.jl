# ---------------------
# elemtype
# 0 -- Node
# 1 —— Line
# 2 —— Triangle
# 3 —— Quadrilateral
# 4 -- Tetrahedron
# 5 -- Hexahedron
# ---------------------
const ET_NUMBER = Dict(
                        Line=>1, 
                        Triangle=>2, 
                        Quadrilateral=>3,
                        Tetrahedron=>4,
                        Hexahedron=>5)
# boundary element type of a specific element type
const BOUND_ET_NUMBER = Dict( 
                            Triangle=>1, 
                            Quadrilateral=>1,
                            Tetrahedron=>2,
                            Hexahedron=>3) 


"""
Supported formats: .msh
"""
function read_gmsh(elemtype::Type{Element{dim,N,M,L}}, msh_file::String) where {dim,N,M,L}

    py_nodeTags, py_nodeCoords = get_nodes(msh_file)
    nnp = length(py_nodeTags)
    # correct the numeration to Julia-style
    nodeTags_order = sortperm(py_nodeTags)
    nodeTags_map = Dict()
    for i = 1:nnp
        nodeTags_map[py_nodeTags[i]] = nodeTags_order[i]
    end
    nodeCoords = py_nodeCoords[1:dim,nodeTags_order]
    reordered_nodeTags = py_nodeTags[nodeTags_order] 
    nodeTags = [k for k = 1:nnp]

    println(py_nodeTags)
    println(nodeCoords)

    # elements -------------------
    py_elemTags, py_elemNodeTags = get_elems(msh_file, ET_NUMBER[elemtype])
    nel = length(py_elemTags)
    for i in eachindex(py_elemNodeTags)
        py_elemNodeTags[i] = nodeTags_map[py_elemNodeTags[i]]
    end
    elemNodeTags = py_elemNodeTags
    # create nodes
    nodes = tags_coords_to_nodes(nodeTags, nodeCoords)

    # elements
    if nel == 0
        error("No element detected!")
    end

    elements = [Element(elemtype, elemNodeTags[1:N,ie]) for ie in 1:nel]
    
    
    # boundary  ---------------------------
    boundElemTags, boundElemNodeTags = get_bounds(msh_file, BOUND_ET_NUMBER[elemtype], dim)

    surface_x = ntuple(m -> ntuple(k -> nodes[boundElemNodeTags[k,m]].x, dim), M)
    surface_u = ntuple(m -> ntuple(k -> nodes[boundElemNodeTags[k,m]].u, dim), M)
    surface_normals = zeros(Float64,dim,M)
    for m in 1:M
        surface_normals[:,m] = getnormal(surface_x[m])
    end
    surface_start, surface_stop = getstartstop(surface_x)
    surface = Surface{M,dim}(surface_x,surface_u, surface_normals, surface_start, surface_stop)

    return nodes, elements, surface
end

function getstartstop(faces::NTuple{M,NTuple{dim,NTuple{dim,Float64}}}) where M where dim
    start, stop = zeros(Float64,dim), zeros(Float64,dim)
    for m = 1:M
        for k = 1:dim
            for axis = 1:dim
                x = faces[m][k][axis] 
                start[axis] = min(x, start[axis])
                stop[axis] = max(x, stop[axis])
            end
        end
    end
    return Tuple(start), Tuple(stop)
end

function tags_coords_to_nodes(nodeTags, nodeCoords)
    nodes = [Node(id,Vec(nodeCoords[:,i]),Vec(nodeCoords[:,i]),Vec(zeros(Float64,dim)),Vec(zeros(Float64,dim)),Vec(zeros(Float64,dim)),Vec(zeros(Float64,dim))) for id in 1:length(nodeTags)]
    return nodes
end