module FiniteElement
using Reexport

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
    time_step!
    

#     J2Plasticity,
#     vonMises,
#     PlasticStructure,
#     Surface,
#     NewtonRaphsonSolver,
#     NewmarkSolver,
#     ExplicitSolver
    
# export
#     generate_grid,
#     add_bc!,
#     
#     fetch_surface,
#     fetch_data,
#     save_to_vtk,
#     time_step!


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

end
