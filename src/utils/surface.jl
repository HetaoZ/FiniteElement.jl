
# ---------------------------------------------
# Immersed surface

"dim 维空间中的 dim-1 维闭合曲面，包含 M 个 face"
struct Surface{M,dim}
    ids::NTuple{M,NTuple{dim,Int}}
    # 存储 M 个 face 的 dim 个结点的 dim 个坐标和速度分量
    x::NTuple{M, NTuple{dim, NTuple{dim,Float64}}}
    u::NTuple{M, NTuple{dim, NTuple{dim,Float64}}}
    # 存储 M 个 face 的法向量的 dim 个坐标（每列代表一个 normal）
    normals::Matrix{Float64}
    # 存储最大和最小坐标
    start::NTuple{dim,Float64} 
    stop::NTuple{dim,Float64}
end

"为了方便单独分析一个 face"
struct SurfaceFace{dim}
    ids::NTuple{dim,Int}
    x::NTuple{dim, NTuple{dim,Float64}}
    u::NTuple{dim, NTuple{dim,Float64}}
    normal::Vector{Float64}
end

"i 是 surface 上的 face 编号"
@inline function getface(surface::Surface{M,dim}, i::Int) where M where dim
    return SurfaceFace{dim}(surface.ids[i], surface.x[i], surface.u[i], surface.normals[:,i])
end



"""
适配 pinpoly 算法
"""
function pinpoly(point::NTuple{dim,Float64}, surface::Surface{M,dim}) where M where dim
    return PointInPoly.pinpoly(point, surface.x, surface.start, surface.stop)
end

# --------------------------------------------------------------
# 针对不同的 package 写接口函数
# --------------------------------------------------------------


# ---------------------------
# for FiniteElement.jl

getsurface!(s::Structure) = getsurface!(s.grid, s.grid.surface_topology)
getsurface!(s::Structure, topo::SurfaceTopology) = getsurface!(s.grid, topo)

function getsurface!(grid::Grid{dim,T}, faces::Matrix{Int}) where {dim,T}
    L, M = size(faces)
    ids = [ntuple(i->grid.nodes[faces[i,i_face]].id, L) for i_face in 1:M]
    x = [ntuple(i->vec2tuple(grid.nodes[faces[i,i_face]].x), L) for i_face in 1:M]
    u = [ntuple(i->vec2tuple(grid.nodes[faces[i,i_face]].u), L) for i_face in 1:M]
    normals = zeros(Float64,dim,size(faces,2))
    for j in axes(normals,2)
        normals[:,j] = getnormal(x[j])
    end
    x = Tuple(x)
    u = Tuple(u)
    ids = Tuple(ids)
    start, stop = getstartstop(x)
    return Surface{M,dim}(ids,x,u,normals,start,stop)
end

function getsurface!(grid::Grid{dim,T}, topo::SurfaceTopology) where {dim,T<:Union{Line,Triangle,Quadrilateral,Tetrahedron}}
    return getsurface!(grid, topo.faces)
end

const ZERO_TRIANGLE = ((0.,0.,0.),(0.,0.,0.),(0.,0.,0.))

function getsurface!(grid::Grid{3,Hexahedron}, topo::SurfaceTopology)
    tri_faces = triangularize(topo.faces)
    return getsurface!(grid, tri_faces)
end

function triangularize(faces::Matrix{Int})
    return hcat(faces[[1,2,3],:], faces[[3,4,1],:])
end

function getquadcenter(nodes::NTuple{4,Node{3}})
    xc = @. (nodes[1].x + nodes[2].x + nodes[3].x + nodes[4].x) / 4
    uc = @. (nodes[1].u + nodes[2].u + nodes[3].u + nodes[4].u) / 4
    return (x = xc, u = uc)
end

function refine_quad(nodes::NTuple{4,Node{3}}, center)
    x = [Tuple.((nodes[1].x, nodes[2].x, center.x)),
         Tuple.((nodes[2].x, nodes[3].x, center.x)),
         Tuple.((nodes[3].x, nodes[4].x, center.x)),
         Tuple.((nodes[4].x, nodes[1].x, center.x))]
    u = [Tuple.((nodes[1].u, nodes[2].u, center.u)),
         Tuple.((nodes[2].u, nodes[3].u, center.u)),
         Tuple.((nodes[3].u, nodes[4].u, center.u)),
         Tuple.((nodes[4].u, nodes[1].u, center.u))]
    start = ntuple(i->minimum((nodes[1].x[i], nodes[2].x[i], nodes[3].x[i], nodes[4].x[i], center.x[i])), 3)
    stop = ntuple(i->maximum((nodes[1].x[i], nodes[2].x[i], nodes[3].x[i], nodes[4].x[i], center.x[i])), 3)
    return x, u, start, stop
end

function getrefinedfaces(quad_nodes, start, stop)
    center = getquadcenter(quad_nodes)
    faces_x, faces_u, faces_start, faces_stop = refine_quad(quad_nodes, center)
    return faces_x, faces_u, ntuple(q->min(start[q], faces_start[q]), 3), ntuple(q->max(stop[q], faces_stop[q]), 3)
end

function getnormal(face::NTuple{2,NTuple{2,Float64}})
    v = [face[2][1] - face[1][1], face[2][2] - face[1][2]]
    n = normalize([v[2], -v[1]])
    return n
end

function getnormal(face::NTuple{3,NTuple{3,Float64}})
    v1 = [face[2][1] - face[1][1], face[2][2] - face[1][2], face[2][3] - face[1][3]]
    v2 = [face[3][1] - face[1][1], face[3][2] - face[1][2], face[3][3] - face[1][3]]
    n = normalize( cross(v1, v2) )
    return n
end

"从 element 的指定 face 提取 surface 所需数据"
function convert_for_surface(face::NTuple{L,Int}, grid::Grid{dim,T}) where {L,dim,T <: Union{Line,Triangle,Quadrilateral,Tetrahedron}}
    x = ntuple(i -> vec2tuple(grid.nodes[face[i]].x), dim)
    u = ntuple(i -> vec2tuple(grid.nodes[face[i]].u), dim)
    start, stop = getstartstop((x,))
    normals = zeros(Float64, dim, 1)
    normals[:,1] = getnormal(x)
    return [x], [u], normals, start, stop
end

"从 element 的指定 face 提取 surface 所需数据"
function convert_for_surface(face::NTuple{L,Int}, grid::Grid{dim,Hexahedron}) where {L,dim}
    start, stop = ntuple(k->0., dim), ntuple(k->0., dim)
    x, u, start, stop = getrefinedfaces([grid.nodes[i] for i in face], start, stop)
    normals = zeros(Float64, dim, length(x))
    for j in eachindex(x)
        normals[:,j] = getnormal(x[j])
    end
    return x, u, normals, start, stop
end

# -------------------------------------
# get surface topology

const EPS = 1e-10
const INF = 1e10

function get_surface_topo!(grid::Grid{1,T}) where T <: AbstractElementType
    if typeof(grid.prototype) <: RectangularGrid  # 方形网格
        x = grid.prototype.start[1]
        start = (x-EPS,)
        stop  = (x+EPS,)
        select_surface!(grid, start, stop)

        x = grid.prototype.stop[1]
        start = (x-EPS,)
        stop  = (x+EPS,)
        select_surface!(grid, start, stop)
    else
        error("Couldn't resolve such grid and get its surface")
    end
end

function get_surface_topo!(grid::Grid{2,T}) where T <: AbstractElementType
    if typeof(grid.prototype) <: RectangularGrid  # 方形网格
        xmin, xmax = grid.prototype.start[1], grid.prototype.stop[1]
        ymin, ymax = grid.prototype.start[2], grid.prototype.stop[2]

        start = (-INF,ymin-EPS)
        stop  = ( INF,ymin+EPS)
        select_surface!(grid, start, stop)

        start = (-INF,ymax-EPS)
        stop  = ( INF,ymax+EPS)
        select_surface!(grid, start, stop)

        start = (xmin-EPS,-INF)
        stop  = (xmin+EPS, INF)
        select_surface!(grid, start, stop)

        start = (xmax-EPS,-INF)
        stop  = (xmax+EPS, INF)
        select_surface!(grid, start, stop)
    else
        error("Couldn't resolve such grid and get its surface")
    end
end

function get_surface_topo!(grid::Grid{3,T}) where T <: AbstractElementType
    if typeof(grid.prototype) <: RectangularGrid  # 方形网格
        xmin, xmax = grid.prototype.start[1], grid.prototype.stop[1]
        ymin, ymax = grid.prototype.start[2], grid.prototype.stop[2]
        zmin, zmax = grid.prototype.start[3], grid.prototype.stop[3]

        start = (-INF,-INF,zmin-EPS)
        stop  = ( INF, INF,zmin+EPS)
        select_surface!(grid, start, stop)

        start = (-INF,-INF,zmax-EPS)
        stop  = ( INF, INF,zmax+EPS)
        select_surface!(grid, start, stop)

        start = (xmin-EPS, -INF, -INF)
        stop  = (xmin+EPS,  INF,  INF)
        select_surface!(grid, start, stop)

        start = (xmax-EPS, -INF, -INF)
        stop  = (xmax+EPS,  INF,  INF)
        select_surface!(grid, start, stop)

        start = (-INF, ymin-EPS, -INF)
        stop  = ( INF, ymin+EPS,  INF)
        select_surface!(grid, start, stop)

        start = (-INF, ymax-EPS, -INF)
        stop  = ( INF, ymax+EPS,  INF)
        select_surface!(grid, start, stop)
    else
        error("Couldn't resolve such grid and get its surface")
    end
end

"将 grid 位于指定 box 空间范围内的 face 加入 surface_topology"
function select_surface!(grid::Grid{dim,T}, start, stop) where {dim,T}

    for elem in grid.elements
        for face in elem.faces
            need_to_add = true
            for node_id in face
                if !point_in_box(grid.nodes[node_id].x, start, stop)
                    need_to_add = false
                end
            end
            if need_to_add
                grid.surface_topology.faces = hcat(grid.surface_topology.faces, reshape(collect(face), (length(face),1)))
            end
        end
    end
end

"将与由 (P1,P2,P3) 三点确定的平面的法向距离小于 bias 的 face 加入 surface_topology"
function select_surface!(grid::Grid{dim,T}, P1, P2, P3; bias = 1e-6) where {dim,T}

    for elem in grid.elements
        for face in elem.faces
            need_to_add = true
            for node_id in face
                if !point_near_plane(grid.nodes[node_id].x, P1, P2, P3, bias = bias)
                    need_to_add = false
                end
            end
            if need_to_add
                grid.surface_topology.faces = hcat(grid.surface_topology.faces, reshape(collect(face), (length(face),1)))
            end
        end
    end
end

"更新 surface 的 normals, start, stop 信息"
function refresh!(surface::Surface{M,dim}) where {M,dim}
    surface.normals = zeros(Float64, dim, M)
    for m in 1:M
        surface.normals[:,m] = getnormal(surface.x[m])
        face_start, face_stop = getstartstop(surface.x[m])
        surface.start = ntuple(i->min(surface.start[i],face_start[i]), dim)
        surface.stop  = ntuple(i->max(surface.stop[i],face_stop[i]), dim)
    end
end

function getstartstop(faces::NTuple{M,NTuple{dim,NTuple{dim,Float64}}}) where M where dim
    start, stop = fill(1e14, dim), fill(-1e14, dim)
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

function point_in_box(point, start, stop)
    return betweeneq(point, start, stop)
end



"不共线的三点确定一个面"
function Plane(P1, P2, P3)
    P1 = collect(P1)
    P2 = collect(P2)
    P3 = collect(P3)
    e1, e2 = P2 - P1, P3 - P1
    normal = normalize(cross(e1, e2))
    d = dot(P1, normal)
    return Plane(normal, d)
end

function point_near_plane(point, plane::Plane; bias = 1e-6)
    @assert bias > 0.
    d = dot(collect(point), plane.normal)
    return abs(d - plane.d) < bias
end

function point_near_plane(point, P1, P2, P3; bias = 1e-6)
    return point_near_plane(point, Plane(P1,P2,P3), bias = bias)
end