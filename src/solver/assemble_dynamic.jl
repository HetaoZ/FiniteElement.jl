"Assembler for dynamic solvers"
function assemble_solution!(solution::TotalLagragianSolution, grid::Grid{dim}, material::AbstractMaterial, states::Vector{Vector{AbstractMaterialState}}, solver::DynamicSolver) where dim
    n_basefuncs = getnbasefunctions(grid.elements[1].ip)
    ndofs = n_basefuncs * dim
    Ke = zeros(Float64, ndofs, ndofs)
    Qe = zeros(Float64, ndofs)
    Me = zeros(Float64, ndofs)
    clear_solution!(solution, solver)

    for (elem, elem_states) in zip(grid.elements, states)
        de = solution.d[getdofs(elem)]
        assemble_element!(elem, grid.nodes, material, de, elem_states, Ke, Qe, Me)
        assemble_global!(solution, elem, Ke, Qe, Me)
    end
end

function clear_solution!(sol::TotalLagragianSolution, ::DynamicSolver)
    fill!(sol.K, 0.)
    fill!(sol.Q, 0.)
    fill!(sol.M, 0.)
end

"""
如果是线弹性，只需要t=0时刻组装一次K和M，后面每次迭代时组装Q。因此这里和静力学求解器有待改进。

load 应该在组装完总体 Q 之后添加到 Q 上，而不是逐个单元加到 Qe 上。

输入参数 Ke, Qe, Me 是为了避免重复分配内存。
"""
function assemble_element!(elem::Element, nodes::Vector{Node{dim}}, material::AbstractMaterial, de::Vector{Float64}, states::Vector{AbstractMaterialState}, Ke::Matrix{Float64}, Qe::Vector{Float64}, Me::Vector{Float64}) where dim

    fill!(Ke, 0.)
    fill!(Qe, 0.)
    fill!(Me, 0.)
    
    x = collect(elem_x(elem, nodes)')
    reinit!(elem.cv, elem.qr, x) # 更新当前的detJ和dNdx
    n = getnbasefunctions(elem.ip)
    nq = length(elem.qr.weights)
    ndofs = dim*n
    
    for i_qpoint in 1:nq
        ε = compute_strain(elem.cv, i_qpoint, de)
        σ, D = compute_stress_tangent(ε, material, states[i_qpoint]) # size(D) = (dim,dim,dim,dim)
        dΩ = getdetJdV(elem.cv, i_qpoint)
        ρₑ = states[i_qpoint].ρ

        # 组装 Ke 和 Me
        ∇N = shape_gradient(elem.cv, i_qpoint) # size(∇N) = (dim,n)
        B = cat(ntuple(k -> node_strain_matrix(∇N[:,k]), n)..., dims = (3,))
        # size(B) = (dim, dim, dim*n)
        for i = 1:ndofs
            Bᵀ_i = Tensor{2,dim}(B[:,:,i]) # 由于B_i对称所以略去转置
            Qe[i] -= (Bᵀ_i ⊡ σ) * dΩ
            for j = 1:ndofs
                B_j = Tensor{2,dim}(B[:,:,j])
                Ke[i,j] += Bᵀ_i  ⊡ D ⊡ B_j * dΩ
                Me[i,j] += Bᵀ_i  ⊡ B_j * dΩ * states[i_qpoint].ρ
            end
        end

        # 集中质量矩阵
        mass = states[i_qpoint].ρ * dΩ
        println("ρ = ", states[i_qpoint].ρ)
        for i = 1:ndofs
            Me[i] += mass
        end
    end
    println("sum(Me) = ", sum(Me))

    # Qe 表示体积力和面力对结点载荷的贡献，此处记为 0， 不做计算，在最终的 Q 上一次性添加外载荷 load

    return Ke, Qe, Me
end


function assemble_global!(sol::TotalLagragianSolution, elem::Element, Ke, Qe, Me)
    global_dofs = getdofs(elem)
    for (i, I) in enumerate(global_dofs)
        for (j, J) in enumerate(global_dofs)
            sol.K[I,J] += Ke[i,j]
        end
        sol.Q[I] += Qe[i]
        sol.M[I,I] += Me[i]
    end
end
