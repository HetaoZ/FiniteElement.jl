function LinearElasticity(dim, E, ν, ρ₀; planar_strain::Bool = true)
    δ(i,j) = i == j ? 1.0 : 0.0 # helper function
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)

    # Isymdev(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
    # 在二维空间中对应于平面应变问题
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    Dᵉ = SymmetricTensor{4, dim}(temp)

    return LinearElasticity{dim, Float64,  typeof(Dᵉ)}(E, ν, G, K, Dᵉ, ρ₀)
end

# function elast_matrix_pstrain(E, ν)
#     E0 = E / (1 - ν^2)
#     ν0 = ν/ (1 - ν)
#     return elast_matrix_pstress(E0, ν0)
# end

# function elast_matrix_pstress(E, ν)
#     D0 = E / (1 - ν^2)
#     D = D0 * [1 ν 0; ν 1 0; 0 0 (1-ν)/2]
#     return D
# end

# function elast_matrix_sym(E, ν)
#     D0 = E*(1-ν) / ((1+ν)*(1-2*ν))
#     a = ν / (1-ν)
#     b = (1-2*ν)/(2*(1-ν))
#     D = D0 * 
#     [1 a 0 a; 
#     a 1 0 a; 
#     0 0 b 0; 
#     a a 0 1]
#     return D
# end

function J2Plasticity(dim, E, ν, σ₀, H, ρ₀)
    δ(i,j) = i == j ? 1.0 : 0.0 # helper function
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)

    Isymdev(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    Dᵉ = SymmetricTensor{4, dim}(temp)

    return J2Plasticity{dim, Float64,  SymmetricTensor{4, dim}}(E, ν, G, K, σ₀, H, Dᵉ, ρ₀)
end



function LinearElasticState(dim)
    return LinearElasticState(
                zero(SymmetricTensor{2, dim}),
                zero(SymmetricTensor{2, dim}),
                0.0,
                0.0)
end

function PlasticState(dim)
    return PlasticState(
                zero(SymmetricTensor{2, dim}),
                zero(SymmetricTensor{2, dim}),
                0.0,
                zero(SymmetricTensor{2, dim}),
                zero(SymmetricTensor{2, dim}),
                0.0,
                0.0,
                0.0)
end

function update_state!(state::LinearElasticState)
    state.σ = state.temp_σ
end

function update_state!(state::PlasticState)
    state.ϵᵖ = state.temp_ϵᵖ
    state.σ = state.temp_σ
    state.k = state.temp_k
end

"""
n = number of elements
nq = number of quadrature points per element
"""
function new_states(::LinearElasticity{dim,T,S}, n::Int, nq::Int) where {dim,T,S}
    return fill(fill(LinearElasticState(dim), nq), n)
end

function new_states(::J2Plasticity{dim,T,S}, n::Int, nq::Int) where {dim,T,S}
    return fill(fill(PlasticState(dim), nq), n)
end

function compute_stress_tangent(ϵ::SymmetricTensor, material::LinearElasticity, state::LinearElasticState)
    # elastic loading
    state.temp_σ = material.Dᵉ ⊡ ϵ
    return state.temp_σ, material.Dᵉ
end

function compute_stress_tangent(ϵ::SymmetricTensor, material::J2Plasticity, state::PlasticState)
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

function symmetrize_lower!(K)
    for i in axes(K, 1)
        for j in i + 1:size(K, 1)
            K[i,j] = K[j,i]
        end
    end
end