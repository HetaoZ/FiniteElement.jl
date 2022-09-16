include("../src/FiniteElement.jl")
using .FiniteElement

nodes, elements, surface = read_gmsh(Triangle, "../in/Triangle.msh")
dump(nodes)