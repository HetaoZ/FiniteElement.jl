function doassemble_K!(cellvalues::CellVectorValues{dim},
    facevalues::FaceVectorValues{dim}, K::SparseMatrixCSC, grid::Grid, dh::DofHandler, material::J2Plasticity, d, states, load) where {dim}

    Q = zeros(ndofs(dh))
    assembler = start_assemble(K, Q)

    nu = getnbasefunctions(cellvalues)
    Qe = zeros(nu)     # element residual vector
    Ke = zeros(nu, nu) # element tangent matrix

    for (cell, state) in zip(CellIterator(dh), states)
        fill!(Ke, 0)
        fill!(Qe, 0)
        eldofs = celldofs(cell)
        de = d[eldofs]
        Ke, Qe = assemble_cell!(Ke, Qe, cell, cellvalues, facevalues, grid, material, de, state, load[eldofs])
        assemble!(assembler, eldofs, Qe, Ke)
    end

    return K, Q
end

function doassemble_M!(M::SparseMatrixCSC, cellvalues::CellVectorValues{dim}, dh::DofHandler, states) where {dim}

    n_basefuncs = getnbasefunctions(cellvalues)
    Me = zeros(n_basefuncs, n_basefuncs)

    assembler = start_assemble(M)

    @inbounds for (cell, state) in zip(CellIterator(dh), states)

        fill!(Me, 0)
        reinit!(cellvalues, cell)

        for q_point in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, q_point)
            ρₑ =  state[q_point].ρ

            for i in 1:n_basefuncs
                v  = shape_value(cellvalues, q_point, i)
                for j in 1:n_basefuncs
                    u = shape_value(cellvalues, q_point, j)
                    Me[i, j] += ρₑ * (v ⋅ u) * dΩ
                end
            end
        end

        assemble!(assembler, celldofs(cell), Me)
    end
    return M
end

function update_state_rho!(grid::Grid, d, states, nqp)
    @inbounds for id in 1:getncells(grid)
        V = volume(grid.cells[id], grid.nodes, d)
        for i = 1:nqp 
            states[id][i].ρ = states[id][i].m₀ / V
        end
    end
end

function doassemble_dynamic!(s::PlasticStructure)
    
    s.system.K, s.system.Q = doassemble_K!(s.cellvalues, s.facevalues, s.system.K, s.grid, s.dh, s.material, s.system.d, s.states, s.load)

    nqp = getnquadpoints(s.cellvalues)
    update_state_rho!(s.grid, s.system.d, s.states, nqp)

    s.system.M = doassemble_M!(s.system.M, s.cellvalues, s.dh, s.states) 
end