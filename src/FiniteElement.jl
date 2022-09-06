module FiniteElement
using Reexport

@reexport using Tensors
using WriteVTK

using LinearAlgebra
using SparseArrays
using ReadGmsh

export
    RectangularGrid,
    Quadrilateral,
    LinearElasticity,
    J2Plasticity,
    NewtonRaphsonSolver,
    NewmarkSolver,
    ExplicitSolver,
    Structure,
    Node

export
    find_nodes,
    add_dirichlet!,
    add_node_force!,
    save_vtk
    

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
#     advance!,
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
include("solver/constrain.jl")
include("structure.jl")


# pre/post-process
include("utils/preprocess.jl")
include("utils/postprocess.jl")

end
