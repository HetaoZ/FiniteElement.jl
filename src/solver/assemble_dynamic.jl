"Assembler for dynamic solvers"
function assemble_solution!(solution::TotalLagragianSolution, grid::Grid{dim}, material::AbstractMaterial, states::Vector{Vector{AbstractMaterialState}}, solver::DynamicSolver) where dim
    n_basefuncs = getnbasefunctions(grid.elements[1].ip)
    ndofs = n_basefuncs * dim
    Ke = zeros(Float64, ndofs, ndofs)
    Qe = zeros(Float64, ndofs)
    Me = zeros(Float64, ndofs)
    clear_solution!(solution, solver)

    for (elem, elem_states) in zip(grid.elements, states)
        fill!(Ke, 0.)
        fill!(Qe, 0.)
        fill!(Me, 0.)
        de = solution.d[getdofs(elem)]
        assemble_element!(elem, grid.nodes, material, de, elem_states, Ke, Qe, Me)
        # assemble_global!(solution, Ke, Qe, Me)
    end
end

function clear_solution!(sol::TotalLagragianSolution, ::DynamicSolver)
    fill!(sol.K, 0.)
    fill!(sol.Q, 0.)
    fill!(sol.M, 0.)
end

"""
load 应该在组装完总体 Q 之后添加到 Q 上，而不是逐个单元加到 Qe 上。

输入参数 Ke, Qe, Me 是为了避免重复分配内存。
"""
function assemble_element!(elem::Quadrilateral, nodes::Vector{Node{dim}}, material::AbstractMaterial, de::Vector{Float64}, states::Vector{AbstractMaterialState}, Ke::Matrix{Float64}, Qe::Vector{Float64}, Me::Vector{Float64}) where dim
    
    x = collect(elem_x(elem, nodes)')
    reinit!(elem.cv, elem.qr, x) # 更新当前的detJ和dNdx
    n = getnbasefunctions(elem.ip)
    nq = length(elem.qr.weights)

    # 组装 Ke
    for i_qpoint in 1:nq
        ε = compute_strain(elem.cv, i_qpoint, de)
        σ, D = compute_stress_tangent(ε, material, states[i_qpoint]) # size(D) = (dim,dim,dim,dim)
        detJdV = getdetJdV(elem.cv, i_qpoint)

        ∇N = shape_gradient(elem.cv, i_qpoint) # size(∇N) = (dim,n)
        B = cat(ntuple(k -> node_strain_matrix(∇N[:,k]), n)..., dims = (3,))
        # size(B) = (dim, dim, dim*n)
        for i = 1:dim*n
            Bᵀ_i = Tensor{2,dim}(B[:,:,i]) # 由于B_i对称所以略去转置
            for j = 1:dim*n
                B_j = Tensor{2,dim}(B[:,:,j])
                Ke[i,j] += Bᵀ_i  ⊡ D ⊡ B_j * detJdV
            end
        end
    end

    # Qe 表示体积力和面力对结点载荷的贡献，此处记为 0， 不做计算，在最终的 Q 上一次性添加外载荷 load

    # 组装 Me


    
    return Ke, Qe, Me
end
