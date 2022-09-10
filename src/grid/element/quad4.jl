"4结点二维四边形单元"
function Element(::Type{Quadrilateral}, connection, faces)
    qr, ip = QuadratureRule{2,RefCube}(2), Lagrange{2,RefCube,1}()
    cv = CellScalarValues(qr, ip)
    return Quadrilateral(connection, faces, cv, qr, ip)
end

function init_volume(elem::Quadrilateral, nodes)
    x = elem_x0(elem, nodes)
    return polygon_area(x[1,:], x[2,:])
end

function volume(elem::Quadrilateral, nodes::Vector{Node})
    x = elem_x(elem, nodes)
    return polygon_area(x[1,:], x[2,:])
end
