"2结点一维线段单元"
function Element(::Type{Line}, connection, faces)
    qr, ip = QuadratureRule{1,RefCube}(2), Lagrange{1,RefCube,1}()
    cv = CellScalarValues(qr, ip)
    return Line(connection, faces, cv, qr, ip)
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
