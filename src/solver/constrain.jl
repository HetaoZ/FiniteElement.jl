# ----------------------------------------------
# 一般约束
function clear_constrains!(s::Structure)
    s.constrains = AbstractConstrain[]
end

# ----------------------------------------------
# 位移边界

function clear_bcs!(s::Structure)
    filter!(c->typeof(c) ∉ (NodeDisplacementConstrain,), s.constrains)
end

function add_bc!(s::Structure, bc::NodeDisplacementConstrain; clear::Bool = false)
    if clear 
        clear_bcs!(s)
    end
    push!(s.constrains, bc)
end

"""
node_ids: 受约束的结点编号向量

constrained_dofs: 每个结点受到约束的自由度，例如 constrained_dofs = [2,3] 表示仅约束第2和第3个自由度

disp_func: 结点位移 d 的受约束分量关于时间 t 的函数，返回一个长度与约束自由度数一致的元组，例如：三维空间中，约束第2和第3个自由度为0，那么 disp_func = t -> (0,0)
"""
function add_bc!(s::Structure{dim}, node_ids::Vector{Int}, constrained_dofs::Vector{Int},
    disp_func::Function; clear::Bool = false) where dim

    @assert length(disp_func(0)) == length(constrained_dofs) && length(constrained_dofs) <= dim

    add_bc!(s, NodeDisplacementConstrain{dim}(node_ids, constrained_dofs, disp_func), clear=clear)
end

# ----------------------------------------------
# 载荷

function clear_forces!(s::Structure)
    filter!(c->typeof(c) ∉ (NodeForceConstrain,), s.constrains)
end

function add_force!(s::Structure{dim}, node_force::NodeForceConstrain{dim}; clear::Bool = false) where dim 
    if clear
        clear_forces!(s)
    end
    push!(s.constrains, node_force)
end

function add_force!(s::Structure{dim}, node_forces::Vector{NodeForceConstrain{dim}}; clear::Bool = false) where dim 
    if clear
        clear_forces!(s)
    end
    append!(s.constrains, node_forces)
end

"""
node_ids: 受约束的结点编号向量

force_func: 结点力f关于坐标和时间(x,t)的函数，返回一个与维数一致的元组，例如：三维空间中 (x,t) -> (100,0,0)
"""
function add_force!(s::Structure{dim}, node_ids::Vector{Int}, force_func::Function; clear::Bool = false) where dim 
    
    @assert length(force_func(0,0)) == dim

    add_force!(s, NodeForceConstrain{dim}(node_ids, force_func), clear = clear)
end

#--------------------------------------------------------------------------
const CONSTRAIN_ALPHA = 1.0e14

"""
施加所有约束
"""
function apply_constrains!(s::Structure, t::Float64)
    for c in s.constrains
        apply_constrain!(s.solution, c, s.grid.nodes, t)
    end
end

"""
施加位移约束
"""
function apply_constrain!(sol::TotalLagragianSolution, constrain::NodeDisplacementConstrain{dim}, nodes::Vector{Node{dim}}, t::Float64) where dim

    for node_id in constrain.node_ids
        for (i, elem_dof) in enumerate(constrain.constrained_dofs)

            global_dof = Int(dim * (node_id - 1) + elem_dof)

            sol.K[global_dof,global_dof] *= CONSTRAIN_ALPHA
            sol.M[global_dof,global_dof] *= CONSTRAIN_ALPHA
            sol.Q[global_dof] = sol.K[global_dof,global_dof] * constrain.disp_func(t)[i]
        end
    end
end

function set_disp!(sol::TotalLagragianSolution, constrain::NodeDisplacementConstrain{dim}, dt, t) where dim
    for node_id in constrain.node_ids
        for (i, elem_dof) in enumerate(constrain.constrained_dofs)
            global_dof = Int(dim * (node_id - 1) + elem_dof)
            
            old_d = sol.d[global_dof]
            old_u = sol.u[global_dof]

            new_d = constrain.disp_func(t)[i]
            new_u = (new_d - old_d) / dt
            new_a = (new_u - old_u) / dt
            
            sol.a[global_dof] = new_a
            sol.u[global_dof] = new_u
            sol.d[global_dof] = new_d
        end
    end
end

function set_disp!(s::Structure, dt, t)
    for c in s.constrains
        if typeof(c) == NodeDisplacementConstrain
            set_disp!(s.sol, c, dt, t)
        end
    end
end

"""
施加力约束
"""
function apply_constrain!(sol::TotalLagragianSolution, constrain::NodeForceConstrain{dim}, nodes::Vector{Node{dim}}, t::Float64) where dim 

    for node_id in constrain.node_ids

        global_dofs = Int(dim * (node_id-1))+1:Int(dim * node_id)
        node_x = nodes[node_id].x

        sol.Q[global_dofs] += collect(constrain.force_func(node_x, t))
    end
end