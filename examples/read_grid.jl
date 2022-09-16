include("../src/FiniteElement.jl")
using .FiniteElement

grid = read_grid("/home/hetao/Projects/Notebooks/triangular_prism.msh")
# display(grid); println()

# 手动划定surface区域
select_surface!(grid, (-1e3,-1e3,-1e3), (1e3,1e3,1e-6))
# 扩展surface区域
select_surface!(grid, (-1e3,-1e3,-1e3), (1e3,1e-6,1e3))
# display(grid.surface_topology.faces); println()

surface = FiniteElement.getsurface!(grid, grid.surface_topology)
# display(surface); println()