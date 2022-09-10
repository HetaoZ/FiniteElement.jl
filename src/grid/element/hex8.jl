"8结点三维六面体单元"
function Element(::Type{Hexahedron}, connection, faces)
    qr, ip = QuadratureRule{3,RefCube}(2), Lagrange{3,RefCube,1}()
    cv = CellScalarValues(qr, ip)
    return Quadrilateral(connection, faces, cv, qr, ip)
end

function init_volume(elem::Hexahedron, nodes)
    x = elem_x0(elem, nodes)
    return hexahedron_volume(x[:,1], x[:,2], x[:,3], x[:,4], x[:,5], x[:,6], x[:,7], x[:,8])
end

function volume(elem::Hexahedron, nodes::Vector{Node})
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