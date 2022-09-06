function TotalLagragianSolution(ndofs::Int)
    d  = zeros(Float64, ndofs)
    Δd = zeros(Float64, ndofs)
    u  = zeros(Float64, ndofs)
    a  = zeros(Float64, ndofs)

    K = spzeros(ndofs, ndofs)
    M = spzeros(ndofs, ndofs)
    # C = spzeros(ndofs, ndofs)

    Q = zeros(Float64, ndofs)
    return TotalLagragianSolution(d,Δd,u,a,K,M,Q)
end

new_solution(s::TotalLagrangianSolver, ndofs::Int) = TotalLagragianSolution(ndofs)