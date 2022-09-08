function doassemble!(s::Structure)
    assemble_solution!(s.solution, s.grid, s.material, s.states, s.solver)
end

function assemble_solution!(solution::TotalLagragianSolution, grid::Grid, material::AbstractMaterial, states::Vector{AbstractMaterialState}, solver::StaticSolver)
    n_basefuncs = getnbasefunctions(grid.elements[1].ip)
    Ke = zeros(Float64, n_basefuncs, n_basefuncs)
    Qe = zeros(Float64, n_basefuncs)
    clear_solution!(solution, solver)

    for (elem, state) in zip(grid.elements, states)
        fill!(Ke, 0.)
        fill!(Qe, 0.)
        de = solution.d[getdofs(elem)]
        assemble_element!(elem, grid.nodes, material, de, state, Ke, Qe)
        # assemble_total!(solution, Ke, Qe)
    end
end

# function assemble_solution!(solution::TotalLagragianSolution, grid::Grid, material::AbstractMaterial, states::Vector{AbstractMaterialState}, solver::DynamicSolver)
#     n_basefuncs = getnbasefunctions(grid.elements[1].ip)
#     Ke = zeros(Float64, n_basefuncs, n_basefuncs)
#     Qe = zeros(Float64, n_basefuncs)
#     Me = zeros(Float64, n_basefuncs)
#     clear_solution!(solution, solver)

#     for (elem, state) in zip(grid.elements, states)
#         fill!(Ke, 0.)
#         fill!(Qe, 0.)
#         fill!(Me, 0.)
#         assemble_element!(elem, grid.nodes, material, de, state, Ke, Qe, Me)
#         assemble_total!(solution, Ke, Qe, Me)
#     end
# end

function clear_solution!(sol::TotalLagragianSolution, ::StaticSolver)
    fill!(sol.K, 0.)
    fill!(sol.Q, 0.)
end

function clear_solution!(sol::TotalLagragianSolution, ::DynamicSolver)
    fill!(sol.K, 0.)
    fill!(sol.Q, 0.)
    fill!(sol.M, 0.)
end

"""
load 应该在组装完总体 Q 之后添加到 Q 上，而不是逐个单元加到 Qe 上。

输入参数 Ke, Qe 是为了避免重复分配内存。
"""
function assemble_element!(elem::Quadrilateral, nodes::Vector{Node{dim}}, material::AbstractMaterial, de::Vector{Float64}, state::AbstractMaterialState, Ke::Matrix{Float64}, Qe::Vector{Float64}) where dim
    
    x = collect(elem_x(elem, nodes)')
    reinit!(elem.cv, elem.ip, elem.qr, x) # 更新当前的detJ和dNdx
    n = getnbasefunctions(elem.ip)
    nq = length(elem.qr.weights)

    for i_qpoint in 1:nq
        ε = compute_strain(elem.cv, i_qpoint, de)
        σ, D = compute_stress_tangent(ε, material, state) # size(D) = (dim,dim,dim,dim)
        detJdV = getdetJdV(elem.cv, i_qpoint)

        ∇N = shape_gradient(elem.cv, i_qpoint) # size(∇N) = (dim,n)
        B = cat(ntuple(k -> node_strain_matrix(∇N[:,k]), n)..., dims = (3,))
        # size(B) = (dim, dim, dim*n)
        Ke +=  transpose(B) ⊡ D ⊡ B * detJdV
    end

    # Qe 表示体积力和面力对结点载荷的贡献，此处记为 0， 不做计算，在最终的 Q 上一次性添加外载荷 load
    
    return Ke, Qe
end

import LinearAlgebra.transpose
function transpose(a::Array{Float64})
    n = length(size(a))
    return permutedims(a, [i for i in n:-1:1])
end

function multiply(a::Array{Float64}, b::Array{Float64})
    sa_left, sa_right = size(a)[1:end-1], size(a)[end]
    sb_left, sb_right = size(b)[1], size(b)[2:end]
    @assert sa_right == sb_left
    return reshape( reshape(a, (prod(sa_left), sa_right)) * 
                    reshape(b, (sb_left, prod(sb_right))) , 
                    (sa_left..., sb_right...))
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
