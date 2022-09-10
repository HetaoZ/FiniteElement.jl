"3结点二维三角形单元"
function Element(::Type{Triangle}, connection, faces)
    qr, ip = QuadratureRule{2,RefTetrahedron}(2), Lagrange{2,RefTetrahedron,1}()
    cv = CellScalarValues(qr, ip)
    return Triangle(connection, faces, cv, qr, ip)
end

function init_volume(elem::Triangle, nodes::Vector{Node})
    x = elem_x0(elem, nodes)
    return polygon_area(x[1,:], x[2,:])
end

function volume(elem::Triangle, nodes::Vector{Node})
    x = elem_x(elem, nodes)
    return polygon_area(x[1,:], x[2,:])
end

function get_min_length(elem::Triangle, nodes)
    return sqrt(volume(elem, nodes) * 2)
end