
advance!(s::PlasticStructure) = s.movable ? static_solver!(s, s.solver) : 0

advance!(s::PlasticStructure, dt::Real) = s.movable ? dynamic_solver!(s, dt, s.solver) : 0

"非线性静力学方程需要多次迭代求解，相当于在无穷小的时间内完成一个驰豫过程。"
function static_solver!(s::PlasticStructure, static_solver::StaticSolver)
    if static_solver == NewtonRaphsonSolver
        newton_raphson_solver!(s, static_solver)
    end
end

function newton_raphson_solver!(s::PlasticStructure, nrsolver::StaticSolver) 
    apply!(s.system.d, s.dbcs)

    newton_itr = -1
    while true; newton_itr += 1

        if newton_itr > 8
            error("Reached maximum Newton iterations, aborting")
            break
        end

        doassemble_static!(s)

        apply_zero!(s.system.K, s.system.Q, s.dbcs)

        norm_r = norm(s.system.Q)

        print("Iteration: $newton_itr \tresidual: $(@sprintf("%.8f", norm_r))\n")
        if norm_r < nrsolver.tolerance
            break
        end

        s.system.Δd = Symmetric(s.system.K) \ s.system.Q
        s.system.d -= s.system.Δd
    end

    for cell_states in s.states
        foreach(update_state!, cell_states)
    end
end

"动力学方程的求解是在有限时间步内的瞬态非平衡解，不需要多次迭代来模拟驰豫过程，因为本来就是非平衡的。"
function dynamic_solver!(s::PlasticStructure, dt::Real, dynamic_solver::DynamicSolver)
    
    apply_zero!(s.system.d, s.dbcs)
    doassemble_dynamic!(s)
    apply_zero!(s.system.K, s.system.Q, s.dbcs)
    core_solver!(s, dt, dynamic_solver.δ, dynamic_solver.α)

    for cell_states in s.states
        foreach(update_state!, cell_states)
    end
end


"需要在System里添加质量矩阵M。先组装好再传入general_solver。"
function core_solver!(s::PlasticStructure, dt::Real, δ::Float64, α::Float64)

    d = s.system.d
    u = s.system.u
    a = s.system.a

    K = Symmetric(s.system.K)
    M = Symmetric(s.system.M)
    Q = s.system.Q
    # vM = sum(Symmetric(s.system.M), dims=1)
    # M = zeros(length(vM), length(vM))
    # for i=1:length(vM)
    #     M[i, i] = vM[i]
    # end

    u_bk = copy(u)
    a_bk = copy(a)

    K_eff = M + α*dt^2*K

    Q_eff = Q - K * (d + dt*u + (0.5-α)*dt^2*a)

    if s.parameters["damping"]
        C = s.system.C
        K_eff += δ*dt*C
        Q_eff -= C * (u + (1-δ)*dt*a)
    end
    
    a = K_eff \ Q_eff

    u = u_bk + ((1-δ)*a_bk + δ*a) * dt
    
    d -= u_bk*dt + ((0.5-α)*a_bk + α*a)*dt^2

    s.system.d = d
    s.system.u = u
    s.system.a = a

    # println("K = ", K)
    # println("M = ", M)
    # println("Q = ", Q)
    # println("d = ", d)

end

function time_step!(s::PlasticStructure)
    if s.solver == NewmarkSolver
        dt = Inf
    else
        error("undef dt")
    end
    return dt
end