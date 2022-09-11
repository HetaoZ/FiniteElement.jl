
# 完全拉格朗日格式
include("total_lagrangian.jl")
include("assemble_static.jl")
include("assemble_dynamic.jl")

# 更新拉格朗日格式（待补充）

# 约束
include("constrain.jl")

"静力学求解"
solve!(s::Structure) = s.movable && static_solver!(s, s.solver)

"动力学求解一个时间步长"
solve!(s::Structure, Δt::Real, t::Real) = s.movable && dynamic_solver!(s, Δt, t, s.solver)

function static_solver!(s::Structure, solver::StaticSolver) 
    if solver == NewtonRaphsonSolver
        newton_raphson_solver!(s, solver)
    end
end


function newton_raphson_solver!(s::Structure, nrsolver::StaticSolver) 
    
    newton_itr = 0
    while true; newton_itr += 1

        if newton_itr > 30
            error("Reached maximum Newton iterations, aborting")
        end

        doassemble!(s)
        apply_constrains!(s, 0.)

        s.solution.Δd =  Symmetric(s.solution.K) \ s.solution.Q  # 需要检验是否符号正确
        s.solution.d += s.solution.Δd 

        for cell_states in s.states
            foreach(update_state!, cell_states)
        end

        update_nodes!(s)

        # ---------------
        # 收敛性检查
        norm_r = norm(s.solution.Q)
        print("Iteration: $newton_itr \tresidual: $(@sprintf("%.8f", norm_r))\n")
        save(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/test2d/structure_"*string(1000+newton_itr))
        if norm_r < nrsolver.tolerance
            break
        end
    end
end

"动力学方程的求解是在有限时间步内的瞬态非平衡解，不需要多次迭代来模拟驰豫过程，因为本来就是非平衡的。"
function dynamic_solver!(s::Structure, dt::Real, t::Real, dynamic_solver::DynamicSolver)
    
    doassemble!(s)
    apply_constrains!(s, t)
    core_solver!(s, dt, dynamic_solver.δ, dynamic_solver.α)
    set_disp!(s, dt, t)

    for cell_states in s.states
        foreach(update_state!, cell_states)
    end

    update_nodes!(s)
end


"需要在solution里添加质量矩阵M。先组装好solution再传入core_solver。"
function core_solver!(s::Structure, dt::Real, δ::Float64, α::Float64)

    u_bk = copy(s.solution.u)
    a_bk = copy(s.solution.a)

    K_eff = s.solution.M + α*dt^2 * s.solution.K
    Q_eff = s.solution.Q - s.solution.K * (s.solution.d + dt * s.solution.u + (0.5-α)*dt^2 * s.solution.a)
    
    s.solution.a = K_eff \ Q_eff
    s.solution.u = u_bk + ((1-δ)*a_bk + δ * s.solution.a) * dt
    s.solution.Δd = - (u_bk*dt + ((0.5-α)*a_bk + α * s.solution.a)*dt^2)
    s.solution.d += s.solution.Δd  

end

time_step!(s::Structure) = time_step!(s, s.solver)
time_step!(::Structure, ::StaticSolver) = Inf

function time_step!(s::Structure, ::DynamicSolver)
    minL = 1.0
    minrho = 1.0

    for e in s.grid.elements        
        minL = get_min_length(e, s.grid.nodes)
        minrho = min(minrho, elem_density(e, s.grid.nodes, s.material))
    end

    C = sqrt(s.material.E/minrho)
    dt = pi * minL / C
    
    return dt
end

function update_nodes!(s::Structure{dim}) where dim
    for id in eachindex(s.grid.nodes)

        node_dofs = (id-1)*dim+1:id*dim

        s.grid.nodes[id].f = Vec(s.solution.Q[node_dofs])

        s.grid.nodes[id].d = Vec(s.solution.d[node_dofs])
        s.grid.nodes[id].x = s.grid.nodes[id].d + s.grid.nodes[id].x0
        if typeof(s.solver) == DynamicSolver
            s.grid.nodes[id].u = Vec(s.solution.u[node_dofs])
            s.grid.nodes[id].a = Vec(s.solution.a[node_dofs])
        end
    end
end