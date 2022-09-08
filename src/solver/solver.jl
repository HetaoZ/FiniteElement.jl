
# 完全拉格朗日格式
include("total_lagrangian.jl")
include("assemble_static.jl")
include("assemble_dynamic.jl")

# 更新拉格朗日格式（待补充）

# 约束
include("constrain.jl")

advance!(s::Structure) = s.movable && static_solver!(s, s.solver)

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

        # apply_constrains!(s.solution, s.constrains)

        norm_r = norm(s.solution.Q)

        print("Iteration: $newton_itr \tresidual: $(@sprintf("%.8f", norm_r))\n")
        if norm_r < nrsolver.tolerance
            break
        end

        s.solution.Δd = -Symmetric(s.solution.K) \ s.solution.Q
        s.solution.d += s.solution.Δd
    end

    for i in eachindex(s.states)
        for j in eachindex(s.states[i])
            update_state!(s.states[i][j])
        end
    end
end