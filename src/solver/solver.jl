
# 完全拉格朗日格式
include("total_lagrangian.jl")
include("assemble_static.jl")
include("assemble_dynamic.jl")

# 更新拉格朗日格式（待补充）

# 约束
include("constrain.jl")

"别名"
advance!(s::Structure) = solve!(s)
advance!(s::Structure, Δt::Real, t::Real) = solve!(s, Δt, t)

"静力学求解"
solve!(s::Structure) = s.movable && static_solver!(s, s.solver)

"动力学求解一个时间步长"
solve!(s::Structure, Δt::Real, t::Real) = s.movable && dynamic_solver!(s, Δt, t, s.solver)

function static_solver!(s::Structure, solver::StaticSolver) 
    newton_raphson_solver!(s, solver)
end

function newton_raphson_solver!(s::Structure, nrsolver::StaticSolver) 
    
    newton_itr = 0
    norm_r1 = 1.0
    while true; newton_itr += 1

        if newton_itr > 50
            error("Reached maximum Newton iterations, aborting")
        end

        doassemble!(s)
        apply_constrains!(s, 0., 0.)

        s.solution.Δd =  Symmetric(s.solution.K) \ s.solution.Q  # 需要检验是否符号正确
        s.solution.d += s.solution.Δd 

        print("itr = ", newton_itr,"  ");display(s.solution.Q);println()


        update_states!(s)
        update_nodes!(s)

        # ---------------
        # 收敛性检查
        norm_r = norm(s.solution.Q)
        if newton_itr == 1
            norm_r1 = norm_r
        end
        print("Iteration: $newton_itr \tresidual: $(@sprintf("%.8f", norm_r/norm_r1))\n")

        save(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/beam3d_static/structure_"*string(1000+newton_itr))

        if norm_r/norm_r1 < nrsolver.tolerance
            break
        end
    end
end

"动力学方程的求解是在有限时间步内的瞬态非平衡解，不需要多次迭代来模拟驰豫过程，因为本来就是非平衡的。"
function dynamic_solver!(s::Structure, dt::Real, t::Real, dynamic_solver::DynamicSolver)
    
    doassemble!(s)
    
    

    apply_constrains!(s, t, dt)
    core_solver3!(s, t, dt, dynamic_solver.δ, dynamic_solver.α)
    core_solver!(s, t, dt, dynamic_solver.δ, dynamic_solver.α)
    # set_disp!(s, dt, t)

    # d = copy(s.solution.d)
    # for i in eachindex(d)
    #     if abs(d[i]) < 1e-14
    #         d[i] = 0.
    #     end
    # end
    # display([s.solution.Q d]);println()
    
    update_states!(s)
    update_nodes!(s)
end


# "两步法。需要在solution里添加质量矩阵M。先组装好solution再传入core_solver。"
function core_solver!(s::Structure, t, dt::Real, δ::Float64, α::Float64)
    println("-- 2-step solver--")

    u_bk = copy(s.solution.u)
    a_bk = copy(s.solution.a)

    K_eff = s.solution.M + α*dt^2 * s.solution.K
    Q_eff = s.solution.Q - s.solution.K * (s.solution.d + dt * s.solution.u + (0.5-α)*dt^2 * s.solution.a)
    
    s.solution.a = K_eff \ Q_eff
    s.solution.u = u_bk + ((1-δ)*a_bk + δ * s.solution.a) * dt
    s.solution.Δd = u_bk*dt + ((0.5-α)*a_bk + α * s.solution.a)*dt^2
    s.solution.d += s.solution.Δd 
    
    dof = 28*3
    println("sol_d = ",  s.solution.d[dof])
    println("sol_u = ",  s.solution.u[dof])
    println("sol_a = ",  s.solution.a[dof])
    # println("sol_Δd = ", s.solution.Δd[dof])
end

"三步法，以位移为变量，可以施加位移边界条件。需要在solution里添加质量矩阵M。先组装好solution再传入core_solver。"
function core_solver3!(s::Structure, t, dt::Real, δ::Float64, α::Float64)
    println()
    println("-- 3-step solver--")

    d_bk = copy(s.solution.d)
    u_bk = copy(s.solution.u)
    a_bk = copy(s.solution.a)

    Q = s.solution.Q
    # if t == 0.
    #     Q_minus_dt = Q
    #     Q_minus_2dt = Q
    # else
        Q_minus_dt = s.solution.Q_minus_dt
        Q_minus_2dt = s.solution.Q_minus_2dt
    # end
    M = s.solution.M
    K = s.solution.K
    d = s.solution.d
    u = s.solution.u
    a = s.solution.a
    d_minus_dt = d - dt * u + dt^2/2*a

    # Q_eff_bar = α*Q + (0.5-2*α+δ)*Q_minus_dt + (0.5-α-δ)*Q_minus_2dt
    Q_eff_bar = α*Q
    Q_eff = Q_eff_bar * dt^2 + (2*M - (0.5-2*α+δ)*dt^2*K) * d + (-M-(0.5+α-δ)*dt^2*K)*d_minus_dt 

    K_eff = M + α*dt^2*K

    sol_d = K_eff \ Q_eff
    sol_a = 1/(α*dt^2) * (sol_d - d_bk) - 1/(α*dt)*u_bk - (1/(2*α) - 1)*a_bk
    sol_u = u_bk + dt*(1-δ)*a_bk + δ*dt*sol_a
    sol_Δd = sol_d - d_bk

    s.solution.Q_minus_2dt = copy(s.solution.Q_minus_dt)
    s.solution.Q_minus_dt  = copy(s.solution.Q)

    dof = 28*3
    println("sol_d = ", sol_d[dof])
    println("sol_u = ", sol_u[dof])
    println("sol_a = ", sol_a[dof])
    # println("sol_Δd = ", sol_Δd[dof])
end

time_step(s::Structure) = time_step(s, s.solver)
time_step(::Structure, ::StaticSolver) = Inf

function time_step(s::Structure, ::DynamicSolver)
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