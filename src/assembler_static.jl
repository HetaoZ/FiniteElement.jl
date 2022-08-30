function compute_stress_tangent(ϵ::SymmetricTensor{2,3}, material::J2Plasticity, state::MaterialState)
    # unpack some material parameters
    G = material.G
    K = material.K
    H = material.H

    # We use (•)ᵗ to denote *trial*-values
    σᵗ = material.Dᵉ ⊡ (ϵ - state.ϵᵖ) # trial-stress
    sᵗ = dev(σᵗ)         # deviatoric part of trial-stress
    J₂ = 0.5 * sᵗ ⊡ sᵗ   # second invariant of sᵗ
    σᵗₑ = sqrt(3.0 * J₂)   # effetive trial-stress (von Mises stress)
    σʸ = material.σ₀ + H * state.k # Previous yield limit

    φᵗ  = σᵗₑ - σʸ # Trial-value of the yield surface

    if φᵗ < 0.0 # elastic loading
        state.temp_σ = σᵗ
        return state.temp_σ, material.Dᵉ
    else # plastic loading
        h = H + 3G
        μ =  φᵗ / h   # plastic multiplier

        c1 = 1 - 3G * μ / σᵗₑ
        s = c1 * sᵗ           # updated deviatoric stress
        σ = s + vol(σᵗ)       # updated stress

        # Compute algorithmic tangent stiffness ``D = \frac{\Delta \sigma }{\Delta \epsilon}``
        κ = H * (state.k + μ) # drag stress
        σₑ = material.σ₀ + κ  # updated yield surface

        δ(i, j) = i == j ? 1.0 : 0.0
        Isymdev(i, j, k, l)  = 0.5 * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k)) - 1.0 / 3.0 * δ(i, j) * δ(k, l)
        Q(i, j, k, l) = Isymdev(i, j, k, l) - 3.0 / (2.0 * σₑ^2) * s[i,j] * s[k,l]
        b = (3G * μ / σₑ) / (1.0 + 3G * μ / σₑ)

        Dtemp(i, j, k, l) = -2G * b * Q(i, j, k, l) - 9G^2 / (h * σₑ^2) * s[i,j] * s[k,l]
        D = material.Dᵉ + SymmetricTensor{4,3}(Dtemp)

        # Store outputs in the material state
        Δϵᵖ = 3 / 2 * μ / σₑ * s            # plastic strain
        state.temp_ϵᵖ = state.ϵᵖ + Δϵᵖ  # plastic strain
        state.temp_k = state.k + μ     # hardening variable
        state.temp_σ = σ               # updated stress
        return state.temp_σ, D
    end
end

function doassemble_static!(cellvalues::CellVectorValues{dim},
    facevalues::FaceVectorValues{dim}, K::SparseMatrixCSC, grid::Grid, dh::DofHandler, material::J2Plasticity, d, states, load) where {dim}
    Q = zeros(ndofs(dh))
    assembler = start_assemble(K, Q)
    # println("-- doassemble: 1--")
    # println(assembler)
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

function doassemble_static!(s::PlasticStructure)
    s.system.K, s.system.Q = doassemble_static!(s.cellvalues, s.facevalues, s.system.K, s.grid, s.dh, s.material, s.system.d, s.states, s.load)
end

function assemble_cell!(Ke, Qe, cell, cellvalues, facevalues, grid, material, de, state, load)
    n_basefuncs = getnbasefunctions(cellvalues)
    reinit!(cellvalues, cell)

    for q_point in 1:getnquadpoints(cellvalues)
        # For each integration point, compute stress and material stiffness
        ∇d = function_gradient(cellvalues, q_point, de)
        ϵ = symmetric(∇d) # Total strain
        σ, D = compute_stress_tangent(ϵ, material, state[q_point])

        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:n_basefuncs
            δϵ = symmetric(shape_gradient(cellvalues, q_point, i))

            Qe[i] += (δϵ ⊡ σ) * dΩ # add internal force to residual
            for j in 1:i
                Δϵ = symmetric(shape_gradient(cellvalues, q_point, j))
                Ke[i, j] += δϵ ⊡ D ⊡ Δϵ * dΩ
            end
        end
    end
    symmetrize_lower!(Ke)

    # 搞清楚这里怎么施加结点力
    Qe -= load

    return Ke, Qe
end

function symmetrize_lower!(K)
    for i in 1:size(K, 1)
        for j in i + 1:size(K, 1)
            K[i,j] = K[j,i]
        end
    end
end;

