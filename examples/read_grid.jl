include("../src/FiniteElement.jl")
using .FiniteElement

grid = read_grid("/home/hetao/Projects/Notebooks/triangular_prism.msh")
# display(grid); println()

# 手动划定surface区域
select_surface!(grid, (0,0,0), (0,0,1), (1,-1,0))
println(size(grid.surface_topology.faces))