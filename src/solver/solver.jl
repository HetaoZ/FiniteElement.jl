
# 完全拉格朗日格式
include("total_lagrangian.jl")
include("assemble_static.jl")
include("assemble_dynamic.jl")

# 更新拉格朗日格式（待补充）

# 约束
include("constrain.jl")

advance!(s::Structure) = s.movable && static_solver!(s, s.solver)
advance!(s::Structure, Δt::Real, t::Real) = s.movable && dynamic_solver!(s, Δt, t, s.solver)

function static_solver!(s::Structure, solver::StaticSolver) 
    if solver == NewtonRaphsonSolver
        newton_raphson_solver!(s, solver)
    end
end

function newton_raphson_solver!(s::Structure, nrsolver::StaticSolver) 
    
    newton_itr = -1
    while true; newton_itr += 1

        if newton_itr > 8
            error("Reached maximum Newton iterations, aborting")
            break
        end

        doassemble!(s)

        apply_constrains!(s, 0.)

        norm_r = norm(s.solution.Q)

        print("Iteration: $newton_itr \tresidual: $(@sprintf("%.8f", norm_r))\n")
        if norm_r < nrsolver.tolerance
            break
        end

        s.solution.Δd = - Symmetric(s.solution.K) \ s.solution.Q  # 需要检验是否符号正确
        s.solution.d += s.solution.Δd 
    end

    for i in eachindex(s.states)
        for j in eachindex(s.states[i])
            update_state!(s.states[i][j])
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
        ext_connection = push!(e.connection, e.connection[1])
        
        for i in eachindex(e.connection)
            L = norm( s.grid.nodes[ext_connection[i]].x - s.grid.nodes[ext_connection[i+1]].x )
            minL = min(minL, L)
        end
    
        minrho = min(minrho, elem_density(e, s.grid.nodes, s.material))
    end

    C = sqrt(s.material.E/minrho)
    dt = pi * minL / C
    
    return dt
end