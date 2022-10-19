
"""
4结点三维四面体单元

结点按如下顺序排列：

第一层：

P3

|  

P1--P2

第二层：

P4
"""
function Element(t::Type{Tetrahedron}, connection)
    quad_rule, ip = QuadratureRule{3,RefTetrahedron}(1), Lagrange{3,RefTetrahedron,1}()
    cv = CellScalarValues(quad_rule, ip)
    faces = infer_faces(t, connection)
    return Tetrahedron(connection, faces, cv, quad_rule, ip)
end

"4结点三维四面体单元"
function Element(::Type{Tetrahedron}, connection, faces)
    quad_rule, ip = QuadratureRule{3,RefTetrahedron}(2), Lagrange{3,RefTetrahedron,1}()
    cv = CellScalarValues(quad_rule, ip)
    return Tetrahedron(connection, faces, cv, quad_rule, ip)
end

"从connection和Element推断faces"
function infer_faces(::Type{Tetrahedron}, c)
    faces = (
    (c[1],c[3],c[2]), 
    (c[1],c[2],c[4]), 
    (c[2],c[3],c[4]), 
    (c[3],c[1],c[4]))
    return faces
end

function init_volume(elem::Tetrahedron, nodes)
    x = elem_x0(elem, nodes)
    return tetrahedron_volume(x[:,1], x[:,2], x[:,3], x[:,4])
end

function volume(elem::Tetrahedron, nodes)
    x = elem_x(elem, nodes)
    return tetrahedron_volume(x[:,1], x[:,2], x[:,3], x[:,4])
end

function get_min_length(elem::Tetrahedron, nodes)
    return cbrt(volume(elem, nodes) * 6)
end

# function tet_minL(P1,P2,P3,P4)
#     min(norm.([P1-P2,P2-P3,P3-P1,P1-P4,P2-P4,P3-P4]))
# end