function TotalLagragianSolution(ndofs::Int)
    d  = zeros(Float64, ndofs)
    Δd = zeros(Float64, ndofs)
    u  = zeros(Float64, ndofs)
    a  = zeros(Float64, ndofs)

    K = spzeros(Float64, ndofs, ndofs)
    M = spzeros(Float64, ndofs, ndofs)
    # C = spzeros(ndofs, ndofs)

    Q = spzeros(Float64, ndofs)
    return TotalLagragianSolution(d,Δd,u,a,K,M,Q)
end

new_solution(s::TotalLagrangianSolver, ndofs::Int) = TotalLagragianSolution(ndofs)