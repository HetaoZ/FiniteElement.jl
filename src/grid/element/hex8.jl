"""
8结点三维六面体单元

结点必须按如下顺序排列：

第一层：

P4--P3

|    |

P1--P2

第二层：

P8--P7

|    |

P5--P6
"""
function Element(et::Type{Hexahedron}, connection)
    qr, ip = QuadratureRule{3,RefCube}(2), Lagrange{3,RefCube,1}()
    cv = CellScalarValues(qr, ip)
    faces = infer_faces(et, connection)
    return Hexahedron(connection, faces, cv, qr, ip)
end

function Element(::Type{Hexahedron}, connection, faces)
    qr, ip = QuadratureRule{3,RefCube}(2), Lagrange{3,RefCube,1}()
    cv = CellScalarValues(qr, ip)
    return Hexahedron(connection, faces, cv, qr, ip)
end

"从connection和Element推断faces"
function infer_faces(::Type{Hexahedron}, c)
    faces = (
    (c[1],c[5],c[8],c[4]), 
    (c[2],c[3],c[7],c[6]), 
    (c[1],c[2],c[6],c[5]), 
    (c[3],c[4],c[8],c[7]), 
    (c[1],c[4],c[3],c[2]), 
    (c[5],c[6],c[7],c[8]))
    return faces
end

function init_volume(elem::Hexahedron, nodes)
    x = elem_x0(elem, nodes)
    return hexahedron_volume(x[:,1], x[:,2], x[:,3], x[:,4], x[:,5], x[:,6], x[:,7], x[:,8])
end

function volume(elem::Hexahedron, nodes::Vector{Node{3}})
    x = elem_x(elem, nodes)
    return hexahedron_volume(x[:,1], x[:,2], x[:,3], x[:,4], x[:,5], x[:,6], x[:,7], x[:,8])
end

function get_min_length(elem::Hexahedron, nodes)
    return cbrt(volume(elem, nodes))
end

# function hex_minL(P1,P2,P3,P4,P5,P6,P7,P8)
#     return min(norm.([P2-P1, P4-P2, P3-P4, P1-P3, P6-P5, P8-P6, P7-P8, P5-P7, P1-P5, P2-P6, P3-P7, P4-P8,
#     P1-P8, P2-P7, P3-P5, P4-P6]))
# end