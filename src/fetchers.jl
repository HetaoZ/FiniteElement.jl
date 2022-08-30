
function fetch_node_coordinate(s::PlasticStructure, i::Int, axis::Int)
    return s.grid.nodes[i].x[axis] + s.system.d[getdim(s.grid)*(i-1)+axis]
end

function fetch_coordinates(s::PlasticStructure)
    return ntuple(axis->ntuple(i->fetch_node_coordinate(s, i, axis), length(s.grid.nodes)), getdim(s.grid))
end

function fetch_load(s::PlasticStructure)
    return s.load
end

function fetch_node_speed(s::PlasticStructure, i::Int, axis::Int)
    return s.system.u[getdim(s.grid)*(i-1)+axis]
end

function fetch_speeds(s::PlasticStructure)
    return ntuple(axis->ntuple(i->fetch_node_speed(s, i, axis), length(s.grid.nodes)), getdim(s.grid))
end

function fetch_surface_x(s::PlasticStructure, dim::Int)
    M = length(boundary_faces)
    cells_x = ntuple(i -> map(node_id -> s.grid.nodes[node_id].x, s.grid.cells[J[i]].nodes), M)
    println("cells_x[1]=",cells_x[1])
end

function fetch_surface(s::PlasticStructure)
    # 根据node的link数无法分辨内外node。
    # 可以根据face的link数分辨内外face，当且仅当link数=1时face在表面上。
    boundary_faces = boundary_matrix_to_faces(s.grid.boundary_matrix, s.grid.cells)
    N = getnnodes(s.grid)
    M = length(boundary_faces)
    dim = getdim(s.grid)
    # 计算outer_normals
    I, J, V = findnz(s.grid.boundary_matrix)
    faces_x = map(face -> map(i -> s.grid.nodes[i].x+s.system.d[dim*(i-1)+1:dim*(i-1)+dim], face), boundary_faces)
    normals = map(face_x -> Vector(get_normal(face_x)), faces_x)
    centers = map(face_x -> get_center(face_x), faces_x)

    # println("node = ", s.grid.nodes[1])
    # println("normals = ", normals[1])
    # println()

    sensors = ntuple(i -> centers[i] + normals[i] * 1.e-10, M)
    cells_x = ntuple(i -> map(node_id -> s.grid.nodes[node_id].x, s.grid.cells[J[i]].nodes), M)
    insides = ntuple(i -> pinpoly(cells_x[i], sensors[i]), M)
    normals = ntuple(i -> insides[i] ? -normals[i] : normals[i], M)

    x = fetch_coordinates(s)
    u = fetch_speeds(s)

    return Surface{N,M,dim}(x, u, boundary_faces, normals)
end

function vec2tuple(v)
    return tuple(v[1], v[2], v[3])
end

function cross_product(a, b)
    if (length(a), length(b)) == (3, 3)
        return eltype(a)[a[2]*b[3]-b[2]*a[3], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-b[1]*a[2]] 
    elseif (length(a), length(b)) == (2, 2)
        return eltype(a)[0,0,a[1]*b[2]-b[1]*a[2]]
    else
        error("undef dim")
    end
end

function pinpoly(cell_x, sensor)
    x = ntuple(i -> cell_x[i][1], 4)
    y = ntuple(i -> cell_x[i][2], 4)
    z = ntuple(i -> cell_x[i][3], 4)
    faces = ((1,2,3), (1,2,4), (2,3,4), (3,1,4))
    return pinpoly(x, y, z, faces, (sensor[1], sensor[2], sensor[3])) == 1
end

function get_center(x)
    return (x[1] + x[2] + x[3])/3
end

function get_normal(x)
    v1 = x[2] - x[1]
    v2 = x[3] - x[1]
    v3 = cross_product(v1, v2)
    normal = v3 / norm(v3)
    return Vec(Tuple(normal))
end

const getfaces = Ferrite.faces

function boundary_matrix_to_faces(m::SparseMatrixCSC, cells)
    I, J, V = findnz(m)
    M = length(I)
    faces = NTuple{3,Int}[]
    for i = 1:M
        face = getfaces(cells[J[i]])[I[i]]
        push!(faces, face)
    end
    return Tuple(faces)
end

function fetch_data(s::PlasticStructure, field)
    @assert field in (:d, :Δd, :u, :a, :f)
    if field == :f
        return fetch_load(s)
    else
        return getfield(s.system, field)
    end
end