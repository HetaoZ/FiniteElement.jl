"4结点三维四面体单元"
function Element(::Type{Tetrahedron}, connection, faces)
    qr, ip = QuadratureRule{3,RefTetrahedron}(2), Lagrange{3,RefCube,1}()
    cv = CellScalarValues(qr, ip)
    return Tetrahedron(connection, faces, cv, qr, ip)
end

function init_volume(elem::Tetrahedron, nodes)
    x = elem_x0(elem, nodes)
    return tetrahedron_volume(x[:,1], x[:,2], x[:,3], x[:,4])
end

function volume(elem::Tetrahedron, nodes::Vector{Node})
    x = elem_x(elem, nodes)
    return tetrahedron_volume(x[:,1], x[:,2], x[:,3], x[:,4])
end

function get_min_length(elem::Tetrahedron, nodes)
    return cbrt(volume(elem, nodes) * 6)
end

# function tet_minL(P1,P2,P3,P4)
#     min(norm.([P1-P2,P2-P3,P3-P1,P1-P4,P2-P4,P3-P4]))
# end