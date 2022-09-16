
"""
4结点二维四边形单元

结点必须按如下顺序排列：

P4--P3

|    |

P1--P2
"""
function Element(t::Type{Quadrilateral}, connection)
    qr, ip = QuadratureRule{2,RefCube}(2), Lagrange{2,RefCube,1}()
    cv = CellScalarValues(qr, ip)
    faces = infer_faces(t, connection)
    return Quadrilateral(connection, faces, cv, qr, ip)
end

function Element(::Type{Quadrilateral}, connection, faces)
    qr, ip = QuadratureRule{2,RefCube}(2), Lagrange{2,RefCube,1}()
    cv = CellScalarValues(qr, ip)
    return Quadrilateral(connection, faces, cv, qr, ip)
end


"从connection和Element推断faces"
function infer_faces(::Type{Quadrilateral}, c)
    faces = (
    (c[1],c[2]), 
    (c[2],c[3]), 
    (c[3],c[4]), 
    (c[4],c[1]))
    return faces
end

function init_volume(elem::Quadrilateral, nodes)
    x = elem_x0(elem, nodes)
    return polygon_area(x[1,:], x[2,:])
end

function volume(elem::Quadrilateral, nodes)
    x = elem_x(elem, nodes)
    return polygon_area(x[1,:], x[2,:])
end

function get_min_length(elem::Quadrilateral, nodes)
    return sqrt(volume(elem, nodes))
end
