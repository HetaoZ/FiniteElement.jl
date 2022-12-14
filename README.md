# FiniteElement

A developing package for 1/2/3D elastic-plastic FEM static/dynamic analysis. DO NOT use it in working environment.

## Installation
```
] add https://github.com/HetaoZ/FiniteElement.jl.git
```

## Examples

Static analysis of 2D elastic beam 
```
using FiniteElement

# define a 2D rectangular grid
nel = (4,2)
grid = RectangularGrid{2,Quadrilateral}((0,0), (2,1), nel)

# define a material
E = 1e9
ν = 0.3
ρ₀ = 1
material = LinearElasticity(2,E,ν,ρ₀)

# define a solver for static analysis
solver = StaticSolver(1.0)

# create a structure
s = Structure(material, grid, solver)

# add displacement boundary conditions
node_ids = find_nodes(s, (-1,-1), (1e-6,1e2))
cdofs = [1,2]
add_disp!(s, node_ids, cdofs, t -> (0,0))

# add external loads
node_ids = find_nodes(s, (2-0.1, 0.5), (2+0.1, 1e2))
add_force!(s, node_ids, (x,t) -> (0,-1e7))

# save data
N = 1000
save(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/test2d/structure_"*string(N))

# solve
solve!(s)

# save final data
save(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/test2d/structure_"*string(N+1))
```
