function create_values(interpolation, qr, face_qr)
    # setup quadrature rules
    # qr      = QuadratureRule{3,RefTetrahedron}(2)
    # face_qr = QuadratureRule{2,RefTetrahedron}(3)

    # create geometric interpolation (use the same as for d)
    interpolation_geom = interpolation

    # cell and facevalues for d
    cellvalues_d = CellVectorValues(qr, interpolation, interpolation_geom)
    facevalues_d = FaceVectorValues(face_qr, interpolation, interpolation_geom)

    return cellvalues_d, facevalues_d
end;

function create_dofhandler(grid, interpolation, dim)
    dh = DofHandler(grid)
    push!(dh, :d, dim, interpolation) # add a displacement field with 3 components
    close!(dh)
    return dh
end

function create_bc(dh, dbc)
    dbcs = ConstraintHandler(dh)
    add!(dbcs, dbc)
    close!(dbcs)
    return dbcs
end;

function add_bc!(s::PlasticStructure, dbc)
    add!(s.dbcs, dbc)
    close!(s.dbcs)
end

function create_grid_handlers(grid::Grid, interpolation::Interpolation, qr, face_qr, dim::Int)
    dh = create_dofhandler(grid, interpolation, dim) # JuaFEM helper function
    cellvalues, facevalues = create_values(interpolation, qr, face_qr)
    return dh, cellvalues, facevalues
end