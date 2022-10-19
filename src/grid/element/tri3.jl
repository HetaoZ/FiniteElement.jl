
"""
3结点二维三角形单元

结点按如下顺序排列：

P3

|  

P1--P2
"""
function Element(t::Type{Triangle}, connection)
    quad_rule, ip = QuadratureRule{2,RefTetrahedron}(1), Lagrange{2,RefTetrahedron,1}()
    cv = CellScalarValues(quad_rule, ip)
    faces = infer_faces(t, connection)
    return Triangle(connection, faces, cv, quad_rule, ip)
end

function Element(::Type{Triangle}, connection, faces)
    quad_rule, ip = QuadratureRule{2,RefTetrahedron}(2), Lagrange{2,RefTetrahedron,1}()
    cv = CellScalarValues(quad_rule, ip)
    return Triangle(connection, faces, cv, quad_rule, ip)
end

"从connection和Element推断faces"
function infer_faces(::Type{Triangle}, c)
    faces = (
    (c[1],c[2]), 
    (c[2],c[3]), 
    (c[3],c[1]))
    return faces
end

function init_volume(elem::Triangle, nodes)
    x = elem_x0(elem, nodes)
    return polygon_area(x[1,:], x[2,:])
end

function volume(elem::Triangle, nodes)
    x = elem_x(elem, nodes)
    return polygon_area(x[1,:], x[2,:])
end

function get_min_length(elem::Triangle, nodes)
    return sqrt(volume(elem, nodes) * 2)
end