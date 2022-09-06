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

condition: 结点位移d关于时间t的函数，例如：t -> Vec(0,0)
"""
function add_bc!(s::Structure{dim}, node_ids::Vector{Int},
    condition::Function; clear::Bool = false) where dim
    add_bc!(s, NodeDisplacementConstrain{dim}(node_ids,condition),clear=clear)
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

condition: 结点力f关于时间t的函数，例如：t -> Vec(0,0)
"""
function add_force!(s::Structure{dim}, node_ids::Vector{Int}, condition::Function; clear::Bool = false) where dim 
    add_force!(s, NodeForceConstrain{dim}(node_ids, condition), clear = clear)
end