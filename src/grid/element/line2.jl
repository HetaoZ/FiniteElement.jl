"2结点一维线段单元"
function Element(et::Type{Line}, connection)
    qr, ip = QuadratureRule{1,RefCube}(2), Lagrange{1,RefCube,1}()
    cv = CellScalarValues(qr, ip)
    faces = infer_faces(et, connection)
    return Line(connection, faces, cv, qr, ip)
end

"2结点一维线段单元"
function Element(::Type{Line}, connection, faces)
    qr, ip = QuadratureRule{1,RefCube}(2), Lagrange{1,RefCube,1}()
    cv = CellScalarValues(qr, ip)
    return Line(connection, faces, cv, qr, ip)
end

"从connection和Element推断faces"
function infer_faces(::Type{Line}, c)
    faces = (
    (c[1],), 
    (c[2],))
    return faces
end

function init_volume(elem::Line, nodes::Vector{Node})
    x = elem_x0(elem, nodes)
    return abs(x[1]-x[2])
end

function volume(elem::Line, nodes::Vector{Node})
    x = elem_x(elem, nodes)
    return abs(x[1]-x[2])
end

function get_min_length(elem::Line, nodes)
    return volume(elem, nodes)
end
