function LinearElasticity(dim, E, ν, ρ₀; planar_strain::Bool = true)
    δ(i,j) = i == j ? 1.0 : 0.0 # helper function
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)

    # Isymdev(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
    # 在二维空间中对应于平面应变问题
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))

    @assert dim in (1,2,3)
    if dim == 3
        Dᵉ = SymmetricTensor{4, dim}(temp)
    elseif dim == 2
        Dᵉ = planar_strain ? elast_matrix_pstrain(E, ν) : elast_matrix_pstress(E, ν)
    else dim == 1
        Dᵉ = [E]
    end

    return LinearElasticity{dim, Float64,  typeof(Dᵉ)}(G, K, Dᵉ, ρ₀)
end

function elast_matrix_pstrain(E, ν)
    E0 = E / (1 - ν^2)
    ν0 = ν/ (1 - ν)
    return elast_matrix_pstress(E0, ν0)
end

function elast_matrix_pstress(E, ν)
    D0 = E / (1 - ν^2)
    D = D0 * [1 ν 0; ν 1 0; 0 0 (1-ν)/2]
    return D
end

function elast_matrix_sym(E, ν)
    D0 = E*(1-ν) / ((1+ν)*(1-2*ν))
    a = ν / (1-ν)
    b = (1-2*ν)/(2*(1-ν))
    D = D0 * 
    [1 a 0 a; 
    a 1 0 a; 
    0 0 b 0; 
    a a 0 1]
    return D
end

function J2Plasticity(dim, E, ν, σ₀, H, ρ₀)
    δ(i,j) = i == j ? 1.0 : 0.0 # helper function
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)

    Isymdev(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    Dᵉ = SymmetricTensor{4, dim}(temp)

    return J2Plasticity{dim, Float64,  SymmetricTensor{4, dim}}(G, K, σ₀, H, Dᵉ, ρ₀)
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

