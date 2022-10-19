"Assembler for structure"
doassemble!(s::Structure) = assemble_solution!(s.solution, s.grid, s.material, s.states, s.solver)

"Assembler for static solvers"
function assemble_solution!(solution::TotalLagragianSolution, grid::Grid{dim}, material::AbstractMaterial, states::Vector{Vector{AbstractMaterialState}}, solver::StaticSolver) where dim
    n_basefuncs = getnbasefunctions(grid.elements[1].ip)
    ndofs = n_basefuncs * dim
    Ke = zeros(Float64, ndofs, ndofs)
    Qe = zeros(Float64, ndofs)
    clear_solution!(solution, solver)

    for i in eachindex(states)
        de = solution.d[getdofs(grid.elements[i])]
        assemble_element(grid.elements[i], grid.nodes, material, de, states[i], Ke, Qe)
        assemble_global!(solution, grid.elements[i], Ke, Qe)
    end

end

function clear_solution!(sol::TotalLagragianSolution, ::StaticSolver)
    fill!(sol.K, 0.)
    fill!(sol.Q, 0.)
end

"""
load 应该在组装完总体 Q 之后添加到 Q 上，而不是逐个单元加到 Qe 上。

输入参数 Ke, Qe 是为了避免重复分配内存。
"""
function assemble_element(elem::Element, nodes::Vector{Node{dim}}, material::AbstractMaterial, de::Vector{Float64}, elem_states::Vector{AbstractMaterialState}, Ke::Matrix{Float64}, Qe::Vector{Float64}) where dim

    fill!(Ke, 0.)
    fill!(Qe, 0.)
    
    x = collect(elem_x(elem, nodes)')
    reinit!(elem.cv, elem.quad_rule, x) # 更新当前的detJ和dNdx
    n = getnbasefunctions(elem.ip)
    nq = length(elem.quad_rule.weights)

    for i_qpoint in 1:nq
        ε = compute_strain(elem.cv, i_qpoint, de)
        compute_stress_tangent!(ε, material, elem_states[i_qpoint]) # size(D) = (dim,dim,dim,dim)
        dΩ = getdetJdV(elem.cv, i_qpoint)

        ∇N = shape_gradient(elem.cv, i_qpoint) # size(∇N) = (dim,n)
        B = cat(ntuple(k -> node_strain_matrix(∇N[:,k]), n)..., dims = (3,))
        # size(B) = (dim, dim, dim*n)
        for i = 1:dim*n
            Bᵀ_i = Tensor{2,dim}(B[:,:,i]) # 由于B_i对称所以略去转置
            Qe[i] -= (Bᵀ_i ⊡ σ) * dΩ
            for j = 1:dim*n
                B_j = Tensor{2,dim}(B[:,:,j])
                Ke[i,j] += Bᵀ_i  ⊡ D ⊡ B_j * dΩ
            end
        end
    end 
end

function assemble_global!(sol::TotalLagragianSolution, elem::Element, Ke, Qe)
    global_dofs = getdofs(elem)
    for (i, I) in enumerate(global_dofs)
        for (j, J) in enumerate(global_dofs)
            sol.K[I,J] += Ke[i,j]
        end
        sol.Q[I] += Qe[i]
    end
end

"""
Compute `Bₖ`

Strain matrix `B = cat(B₁,B₂,...,dims=(3,))`
"""
function node_strain_matrix(v::Vector{Float64})
    dim = length(v)
    Bₖ = zeros(Float64,dim,dim,dim)
    vhalf = v * 0.5
    for i = 1:dim
        Bₖ[i,:,:] = diagm(fill(vhalf[i],dim))
        Bₖ[i,:,i] += vhalf
    
        # e.g., when dim = 3
        # Bₖ[1,1,:] = [v[1],   0,   0]
        # Bₖ[1,2,:] = [v[2],v[1],   0]*0.5
        # Bₖ[1,3,:] = [v[3],   0,v[1]]*0.5

        # Bₖ[2,1,:] = [v[2],v[1],   0]*0.5
        # Bₖ[2,2,:] = [0,   v[2],   0]
        # Bₖ[2,3,:] = [0,   v[3],v[2]]*0.5

        # Bₖ[3,1,:] = [v[3],   0,v[1]]*0.5
        # Bₖ[3,2,:] = [0,   v[3],v[2]]*0.5
        # Bₖ[3,3,:] = [0,      0,v[3]]
    end
    return Bₖ
end

function compute_strain(cv::CellValues{dim}, i_qpoint::Int, de::Vector{Float64}) where dim
    d = collect(reshape(de, dim, Int(length(de)/dim))')
    ∇d = disp_gradient(cv, i_qpoint, d)
    ε = convert(SymmetricTensor, 0.5 * (∇d + transpose(∇d)))
    return ε
end
