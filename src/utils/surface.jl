
# ---------------------------------------------
# Immersed surface

"dim 维空间中的 dim-1 维闭合曲面，包含 M 个 face"
struct Surface{M,dim}
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
    x::NTuple{dim, NTuple{dim,Float64}}
    u::NTuple{dim, NTuple{dim,Float64}}
end

"i 是 surface 上的 face 编号"
@inline function getface(surface::Surface{M,dim}, i::Int) where M where dim
    return SurfaceFace{dim}(surface.x[i], surface.u[i])
end

using PointInPoly

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

function getsurface!(s::FiniteElement.Structure{1})
    if typeof(s.grid.prototype) <: FiniteElement.RectangluarGrid  # 方形网格
        x = ((s.grid.nodes[1].x,), (s.grid.nodes[2].x,))
        u = ((s.grid.nodes[1].u,), (s.grid.nodes[2].u,))
        normals = zeros(Float64, 1, 1)
        normals[:,1] = [-1]
        normals[:,2] = [1]
        start, stop = s.grid.nodes[1].x, s.grid.nodes[2].x
        return Surface{M,1}(x, u, normals, start, stop)
    else
        error("Couldn't resolve such grid and get its surface")
    end
end

function getsurface!(s::FiniteElement.Structure{2})
    if typeof(s.grid.prototype) <: FiniteElement.RectangluarGrid  # 方形网格
        nel = s.grid.prototype.nel
        M = sum(nel) * 2  # n faces
        nodes = reshape(s.grid.nodes, nel .+ 1)
        x = fill(((0.,0.), (0.,0.)), M)
        u = fill(((0.,0.), (0.,0.)), M)
        xmin, xmax, ymin, ymax = nodes[1].x[1], nodes[end].x[1], nodes[1].x[2], nodes[end].x[2]
        for i = 1:nel[1]
            # Bottom: nodes[1:nx+1,1]
            x[i] = Tuple.((nodes[i,1].x, nodes[i+1,1].x))
            u[i] = Tuple.((nodes[i,1].u, nodes[i+1,1].u))
            xmin = min(xmin, nodes[i,1].x[1])
            xmax = max(xmax, nodes[i,1].x[1])
            ymin = min(ymin, nodes[i,1].x[2])
            ymax = max(ymax, nodes[i,1].x[2])
            # Top: nodes[1:nx+1,end]
            x[i+nel[1]] = Tuple.((nodes[i,end].x, nodes[i+1,end].x))
            u[i+nel[1]] = Tuple.((nodes[i,end].u, nodes[i+1,end].u))
            xmin = min(xmin, nodes[i,end].x[1])
            xmax = max(xmax, nodes[i,end].x[1])
            ymin = min(ymin, nodes[i,end].x[2])
            ymax = max(ymax, nodes[i,end].x[2])
        end
        for j = 1:nel[2]
            # Left: nodes[1,1:ny+1]
            x[i+2*nel[1]] = Tuple.((nodes[1,j].x, nodes[1,j+1].x))
            u[i+2*nel[1]] = Tuple.((nodes[1,j].u, nodes[1,j+1].u))
            xmin = min(xmin, nodes[1,j].x[1])
            xmax = max(xmax, nodes[1,j].x[1])
            ymin = min(ymin, nodes[1,j].x[2])
            ymax = max(ymax, nodes[1,j].x[2])
            # Right: nodes[end,1:ny+1]
            x[i+2*nel[1]+nel[2]] = Tuple.((nodes[end,j].x, nodes[end,j+1].x))
            u[i+2*nel[1]+nel[2]] = Tuple.((nodes[end,j].u, nodes[end,j+1].u))
            xmin = min(xmin, nodes[end,j].x[1])
            xmax = max(xmax, nodes[end,j].x[1])
            ymin = min(ymin, nodes[end,j].x[2])
            ymax = max(ymax, nodes[end,j].x[2])
        end
        normals = zeros(Float64, 2, M)
        for (i,face) in enumerate(x)
            normals[:,i] = getnormal(face)
        end
        start, stop = (xmin, ymin), (xmax, ymax)
        return Surface{M,2}(Tuple(x), Tuple(u), normals, start, stop)
    else
        error("Couldn't resolve such grid and get its surface")
    end
end

function getsurface!(s::FiniteElement.Structure{3})
    if typeof(s.grid.prototype) <: FiniteElement.RectangularGrid  # 方形网格
        nel = s.grid.prototype.nel
        M = (nel[1]*nel[3] + nel[2]*nel[3] + nel[1]*nel[2])*8 # M faces (refined)
        nodes = reshape(s.grid.nodes, nel .+ 1)
        x = fill(((0.,0.,0.), (0.,0.,0.),(0.,0.,0.)), M)
        u = fill(((0.,0.,0.), (0.,0.,0.),(0.,0.,0.)), M)
        start = nodes[1].x
        stop = nodes[end].x

        for i = 1:nel[1], k = 1:nel[3]
            # South: nodes[:,1,:]
            i_quad = i+(k-1)*nel[1]
            faces_range = i_quad*4-3:i_quad*4

            x[faces_range], u[faces_range], start, stop = getrefinedfaces((nodes[i,1,k], 
                             nodes[i+1,1,k], 
                             nodes[i+1,1,k+1], 
                             nodes[i,1,k+1]), 
                             start, stop)
 
            # North: nodes[:,end,:]
            i_quad = i+(k-1)*nel[1] + nel[1]*nel[3]
            faces_range = i_quad*4-3:i_quad*4

            x[faces_range], u[faces_range], start, stop = getrefinedfaces((nodes[i,end,k], 
                             nodes[i,end,k+1], 
                             nodes[i+1,end,k+1], 
                             nodes[i+1,end,k]), 
                             start, stop)
        end
        for j = 1:nel[2], k = 1:nel[3]
            # West: nodes[1,:,:]
            i_quad = j+(k-1)*nel[2] + 2*nel[1]*nel[3]
            faces_range = i_quad*4-3:i_quad*4

            x[faces_range], u[faces_range], start, stop = getrefinedfaces((nodes[1,j,k], 
                             nodes[1,j,k+1], 
                             nodes[1,j+1,k+1], 
                             nodes[1,j+1,k]), 
                             start, stop)
 
            # East: nodes[end,:,:]
            i_quad = j+(k-1)*nel[2] + 2*nel[1]*nel[3] + nel[2]*nel[3]
            faces_range = i_quad*4-3:i_quad*4

            x[faces_range], u[faces_range], start, stop = getrefinedfaces((nodes[end,j,k], 
                             nodes[end,j+1,k], 
                             nodes[end,j+1,k+1], 
                             nodes[end,j,k+1]), 
                             start, stop)
        end
        for i = 1:nel[1], j = 1:nel[2]
            # Bottom: nodes[:,:,1]
            i_quad = i+(j-1)*nel[1] + 2*(nel[1]+nel[2])*nel[3]
            faces_range = i_quad*4-3:i_quad*4

            x[faces_range], u[faces_range], start, stop = getrefinedfaces((nodes[i,j,1], 
                             nodes[i,j+1,1], 
                             nodes[i+1,j+1,1], 
                             nodes[i+1,j,1]), 
                             start, stop)
 
            # Top: nodes[:,:,end]
            i_quad = i+(j-1)*nel[1] + 2*(nel[1]+nel[2])*nel[3] + nel[1]*nel[2]
            faces_range = i_quad*4-3:i_quad*4

            x[faces_range], u[faces_range], start, stop = getrefinedfaces((nodes[i,j,end], 
                             nodes[i+1,j,end], 
                             nodes[i+1,j+1,end], 
                             nodes[i,j+1,end]), 
                             start, stop)
        end
        normals = zeros(Float64, 3, M)
        for (i,face) in enumerate(x)
            normals[:,i] = getnormal(face)
        end
        return Surface{M,3}(Tuple(x), Tuple(u), normals, start, stop)
    else
        error("Couldn't resolve such grid and get its surface")
    end
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