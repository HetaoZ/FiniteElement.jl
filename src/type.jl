
## ----------------------------------------------------------
# 材料
# ----------------------------------------------------------
abstract type AbstractMaterial end

struct LinearElasticity{dim,T,S} <: AbstractMaterial
    G::T  # Shear modulus
    K::T  # Bulk modulus
    Dᵉ::S # Elastic stiffness tensor
    ρ₀::T # Initial density
end

struct J2Plasticity{dim, T, S} <: AbstractMaterial
    G::T  # Shear modulus
    K::T  # Bulk modulus
    σ₀::T # Initial yield limit
    H::T  # Hardening modulus
    Dᵉ::S # Elastic stiffness tensor
    ρ₀::T # Initial density
end

abstract type AbstractMaterialState{T,S} end

mutable struct LinearElasticState{T,S} <: AbstractMaterialState{T,S}
    σ::S # stress
    temp_σ::S

    m₀::T
    ρ::T
end

"材料非线性，包括塑性、非线性弹性和粘性"
mutable struct PlasticState{T,S} <: AbstractMaterialState{T,S}
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

# ----------------------------------------------------------
# 网格
# ----------------------------------------------------------

mutable struct Node{dim}
    id::Int
    x0::Vec{dim,Float64}
    x::Vec{dim,Float64}
    d::Vec{dim,Float64}
    u::Vec{dim,Float64}
    a::Vec{dim,Float64}
    f::Vec{dim,Float64}
end

abstract type AbstractElementType{dim,N,M,L} end

"""
Node数据集中存储在Grid里，而单元只存储拓扑信息。

dim: dimension

N: number of element nodes

M: number of element faces

L: number of nodes in each face
"""
struct Element{dim,N,M,L} <: AbstractElementType{dim,N,M,L}
    connection::NTuple{N,Int}
    faces::NTuple{M,NTuple{L,Int}}
end

# 细分的单元类型
const Line = Element{1,2,2,1}
const Triangle = Element{2,3,3,2}
const Quadrilateral = Element{2,4,4,2}
const Tetrahedron = Element{3,4,4,3}
const Hexahedron = Element{3,8,6,4}

"""
Represents a reference shape which quadrature rules and interpolations are defined on.
Currently, the only concrete types that subtype this type are `RefCube` in 1, 2 and 3 dimensions,
and `RefTetrahedron` in 2 and 3 dimensions.
"""
abstract type AbstractRefShape end

struct RefTetrahedron{dim} <: AbstractRefShape end
struct RefCube{dim} <: AbstractRefShape end

struct QuadratureRule{dim,shape<:AbstractRefShape,T}
    weights::Vector{T}
    points::Vector{Vec{dim,T}}
end

"""
    Interpolation{dim, ref_shape, order}()

Return an `Interpolation` of given dimension `dim`, reference shape
(see see [`AbstractRefShape`](@ref)) `ref_shape` and order `order`.
`order` corresponds to the highest order term in the polynomial.
The interpolation is used to define shape functions to interpolate
a function between nodes.

The following interpolations are implemented:

* `Lagrange{1,RefCube,1}`
* `Lagrange{1,RefCube,2}`
* `Lagrange{2,RefCube,1}`
* `Lagrange{2,RefCube,2}`
* `Lagrange{2,RefTetrahedron,1}`
* `Lagrange{2,RefTetrahedron,2}`
* `Lagrange{3,RefCube,1}`
* `Serendipity{2,RefCube,2}`
* `Serendipity{3,RefCube,2}`
* `Lagrange{3,RefTetrahedron,1}`
* `Lagrange{3,RefTetrahedron,2}`

# Examples
```jldoctest
julia> ip = Lagrange{2,RefTetrahedron,2}()
Ferrite.Lagrange{2,Ferrite.RefTetrahedron,2}()

julia> getnbasefunctions(ip)
6
```
"""
abstract type Interpolation{dim,shape,order} end

struct Lagrange{dim,shape,order} <: Interpolation{dim,shape,order} end

"""
Abstract type which has `CellValues` and `FaceValues` as subtypes
"""
abstract type Values{dim,T,refshape} end
abstract type CellValues{dim,T,refshape} <: Values{dim,T,refshape} end
abstract type FaceValues{dim,T,refshape} <: Values{dim,T,refshape} end

struct CellScalarValues{dim,T<:Real,refshape<:AbstractRefShape} <: CellValues{dim,T,refshape}
    N::Matrix{T}
    dNdx::Matrix{Vec{dim,T}}
    dNdξ::Matrix{Vec{dim,T}}
    detJdV::Vector{T}
    M::Matrix{T}
    dMdξ::Matrix{Vec{dim,T}}
    qr::QuadratureRule{dim,refshape,T}
    # The following fields are deliberately abstract -- they are never used in
    # performance critical code, just stored here for convenience.
    func_interp::Interpolation{dim,refshape}
    geo_interp::Interpolation{dim,refshape}
end



abstract type AbstractGrid{dim} end

mutable struct Grid{dim, T<:AbstractElementType} <: AbstractGrid{dim}
    nodes::Vector{Node{dim}}
    elements::Vector{T}
end

abstract type AbstractGridPrototype end

struct RectangularGrid{dim,T<:AbstractElementType} <: AbstractGridPrototype
    start
    stop
    nel::NTuple{dim,Int}
end

# ----------------------------------------------------------
# 求解格式
# ----------------------------------------------------------

abstract type AbstractSolver end
abstract type AbstractSolution end

abstract type TotalLagrangianSolver <: AbstractSolver end

struct StaticSolver <: TotalLagrangianSolver 
    tolerance::Float64
end

const NewtonRaphsonSolver = StaticSolver(1.0)

struct DynamicSolver <: TotalLagrangianSolver 
    δ::Float64
    α::Float64
end

const NewmarkSolver = DynamicSolver(0.52, 0.25*(0.5 + 0.52)^2) # δ = 0.52, α = 0.25*(0.5+δ)^2
const ExplicitSolver = DynamicSolver(0.5, 0.0)

mutable struct TotalLagragianSolution <: AbstractSolution
    # 能否把矢量也写成稀疏形式？
    d::Vector{Float64}
    Δd::Vector{Float64}
    u::Vector{Float64}
    a::Vector{Float64}

    # mises_values::Vector{Float64}
    # κ_values::Vector{Float64}

    K::SparseMatrixCSC
    M::SparseMatrixCSC
    # C::SparseMatrixCSC  # damping matrix 阻尼矩阵

    Q::Vector{Float64}
end

"""
位移和应力边界条件均可表示为结点约束
"""
abstract type AbstractConstrain{dim} end

"""
Dirichlet条件，输入 t，返回一个约束后的变量

施加在点、面、体上的力载荷均可以视为应力边界条件，进而等效为结点的Dirichlet约束，例如：

condition = (node::Node{3}, t) -> begin
    node.f = Vec(1e3,1e3,0)
    node 
end
"""
abstract type NodeConstrain{dim} <: AbstractConstrain{dim} end

"""
位移边界条件

dim: 维数

node_ids: 受约束的结点编号向量

condition: 结点位移d关于时间t的函数，例如：t -> Vec(0,0)
"""
struct NodeDisplacementConstrain{dim} <: NodeConstrain{dim}
    node_ids::Vector{Int}
    condition::Function
end

"""
力边界条件

dim: 维数

node_ids: 受约束的结点编号向量

condition: 结点力f关于时间t的函数，例如：t -> Vec(0,0)
"""
struct NodeForceConstrain{dim} <: NodeConstrain{dim}
    node_ids::Vector{Int}
    condition::Function
end

# ----------------------------------------------
# 结构
# ----------------------------------------------

mutable struct Structure{dim}
    material::AbstractMaterial
    grid::AbstractGrid
    solver::AbstractSolver
    solution::AbstractSolution
    constrains::Vector{AbstractConstrain}
    parameters::Dict
    movable::Bool
end