include("../src/FiniteElement.jl")
using .FiniteElement

grid = read_grid("/home/hetao/Projects/Notebooks/triangular_prism.msh")

# display(grid); println()

surface = extend_surface(grid, (-1e3,-1e3,-1e3), (1e3,1e3,1e-6))

display(surface); println()