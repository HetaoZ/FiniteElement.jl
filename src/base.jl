struct SurfaceNode{dim}
    id::Int
    x::Vector{Float64}
    u::Vector{Float64}
    function SurfaceNode(id::Int, x::Vector{Float64}, u::Vector{Float64})
        new{length(x)}(id, x, u)
    end
end

struct SurfaceFace{dim}
    id::Int
    nodes::NTuple{dim,SurfaceNode}
    normal::Vector{Float64}
end

function getarea(face::SurfaceFace{1})
    return 1.0
end

function getarea(face::SurfaceFace{2})
    return norm(face.nodes[1].x - face.nodes[2].x)
end

function getarea(face::SurfaceFace{3})
    A, B, C = face.nodes[1].x, face.nodes[2].x, face.nodes[3].x
    a, b, c = norm(A - B), norm(B - C), norm(C - A)
    p = (a+b+c)/2
    return sqrt(p*(p-a)*(p-b)*(p-c))
end

"boundary polygon/polyhedron surface of N nodes, M faces, dim dimension"
struct Surface{N,M,dim}
    x::NTuple{dim, NTuple{N,Float64}}
    u::NTuple{dim, NTuple{N,Float64}}
    faces::NTuple{M, NTuple{dim,Int}} 
    # 只需要判定探针位于当前cell内外即可，因为已知这些faces均为表面face，所以探针位于当前cell外是探针位于表面外的充分必要条件。
    normals::NTuple{M, Vector{Float64}}
end

@inline function getnode(surface::Surface, i::Int)
    return SurfaceNode(i, [coord[i] for coord in surface.x], [speed[i] for speed in surface.u])
end

@inline function getface(surface::Surface{N,M,dim}, i::Int) where N where M where dim
    faceids = surface.faces[i]
    return SurfaceFace(i, ntuple(k -> getnode(surface, faceids[k]), dim), Vector(surface.normals[i]))
end

abstract type AbstractPlasticity end

struct J2Plasticity{T, S <: SymmetricTensor{4, 3, T}} <: AbstractPlasticity
    G::T  # Shear modulus
    K::T  # Bulk modulus
    σ₀::T # Initial yield limit
    H::T  # Hardening modulus
    Dᵉ::S # Elastic stiffness tensor
    ρ₀::T # Initial density
end

function J2Plasticity(E, ν, σ₀, H, ρ₀)
    δ(i,j) = i == j ? 1.0 : 0.0 # helper function
    G = E / 2(1 + ν)
    K = E / 3(1 - 2ν)

    Isymdev(i,j,k,l) = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
    temp(i,j,k,l) = 2.0G *( 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) + ν/(1.0-2.0ν)*δ(i,j)*δ(k,l))
    Dᵉ = SymmetricTensor{4, 3}(temp)
    return J2Plasticity(G, K, σ₀, H, Dᵉ, ρ₀)
end

mutable struct MaterialState{T, S <: SecondOrderTensor{3, T}}
    # Store "converged" values
    ϵᵖ::S # plastic strain
    σ::S # stress
    k::T # hardening variable

    # Store temporary values used during equilibrium iterations
    temp_ϵᵖ::S
    temp_σ::S
    temp_k::T

    m₀::T
    ρ::T
end

function MaterialState()
    return MaterialState(
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                0.0,
                zero(SymmetricTensor{2, 3}),
                zero(SymmetricTensor{2, 3}),
                0.0,
                0.0,
                0.0)
end

function update_state!(state::MaterialState)
    state.ϵᵖ = state.temp_ϵᵖ
    state.σ = state.temp_σ
    state.k = state.temp_k
end

function vonMises(σ, dim)
    if dim == 3
        s = dev(σ)
        return sqrt(3.0/2.0 * s ⊡ s)
    else
        error("undef")
    end
end

abstract type AbstractSolver end

struct StaticSolver <: AbstractSolver 
    tolerance::Float64
end

const NewtonRaphsonSolver = StaticSolver(1.0)

struct DynamicSolver <: AbstractSolver 
    δ::Float64
    α::Float64
end

const NewmarkSolver = DynamicSolver(0.52, 0.25*(0.5 + 0.52)^2) # δ = 0.52, α = 0.25*(0.5+δ)^2
const ExplicitSolver = DynamicSolver(0.5, 0.0)

mutable struct PlasticSystem
    # 能否把矢量也写成稀疏形式？
    d::Vector{Float64}
    Δd::Vector{Float64}
    u::Vector{Float64}
    a::Vector{Float64}

    # mises_values::Vector{Float64}
    # κ_values::Vector{Float64}

    K::SparseMatrixCSC
    M::SparseMatrixCSC
    C::SparseMatrixCSC

    Q::Vector{Float64}
end

function PlasticSystem(dh::DofHandler)
    n_dofs = ndofs(dh)
    d  = zeros(n_dofs)
    Δd = zeros(n_dofs)
    u  = zeros(n_dofs)
    a  = zeros(n_dofs)
    
    K = create_sparsity_pattern(dh)
    M = create_sparsity_pattern(dh)
    C = create_sparsity_pattern(dh)

    Q  = zeros(n_dofs)

    return PlasticSystem(d, Δd, u, a, K, M, C, Q)
end

mutable struct PlasticStructure
    dim::Int
    material::AbstractPlasticity
    grid::Grid
    dh::DofHandler
    dbcs::ConstraintHandler
    cellvalues::CellVectorValues
    facevalues::FaceVectorValues
    states::Vector{Vector{MaterialState}}
    system::PlasticSystem
    load::Vector{Float64}
    parameters::Dict
    solver::AbstractSolver
    movable::Bool
    # constrained_dofs::Vector{Int}
end

function PlasticStructure(material::AbstractPlasticity, grid::Grid, interpolation, qr, face_qr, dim::Int)

    dh, cellvalues, facevalues = create_grid_handlers(grid, interpolation, qr, face_qr, dim) 
    
    dbcs = ConstraintHandler(dh)

    nqp = getnquadpoints(cellvalues)
    states = [[MaterialState() for _ in 1:nqp] for _ in 1:getncells(grid)]

    system = PlasticSystem(dh)

    n_dofs = ndofs(dh)
    load = zeros(n_dofs)

    for id in 1:getncells(grid)
        V₀ = volume(grid.cells[id], grid.nodes, system.d)
        m₀ = V₀ * material.ρ₀
        for i = 1:nqp 
            states[id][i].m₀ = m₀
            states[id][i].ρ = material.ρ₀
        end
    end

    parameters = Dict(
        "ndofs"=> n_dofs,
        "nnodes"=>getnnodes(grid),
        "damping"=>false
    )
    
    return PlasticStructure(dim, material, grid, dh, dbcs, cellvalues, facevalues, states, system, load, parameters, NewmarkSolver, true)
end

# PlasticStructure(material, grid, interpolation, qr, face_qr) = PlasticStructure(material, grid, interpolation, qr, face_qr, interpolation)

# function Base.copy(s::PlasticStructure)

# end
