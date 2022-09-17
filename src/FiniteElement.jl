module FiniteElement
using Reexport

include("/home/hetao/Projects/JuliaProj/PointInPoly.jl/src/PointInPoly.jl")
using .PointInPoly
using Tensors
using WriteVTK
using LinearAlgebra
using SparseArrays
using Base: @propagate_inbounds
using Printf

export
    RectangularGrid,
    Line,
    Triangle,
    Quadrilateral,
    Tetrahedron,
    Hexahedron,
    LinearElasticity,
    J2Plasticity,
    NewtonRaphsonSolver,
    StaticSolver,
    NewmarkSolver,
    ExplicitSolver,
    DynamicSolver,
    Structure

export
    solve!,
    advance!,
    find_nodes,
    add_bc!,
    add_force!,
    save,
    time_step!,
    read_grid,
    select_surface!
    


# main files
include("type.jl")
include("utils/utils.jl")
include("material/material.jl")
include("grid/grid.jl")
include("solver/solver.jl")

# pre/post-process
include("utils/preprocess.jl")
include("utils/postprocess.jl")
include("utils/surface.jl")
include("utils/FiniteMesh_interface.jl")

end
