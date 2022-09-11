include("../src/FiniteElement.jl")
using .FiniteElement

# define a 2D rectangular grid
nel = (2,1,1)
grid = RectangularGrid{3, Hexahedron}((0,0,0), (2,1,1), nel)


# define a material
E = 1e9
ν = 0.3
ρ₀ = 1
material = LinearElasticity(3,E,ν,ρ₀)

# define a solver for static analysis
solver = StaticSolver(1e-3)

# create a structure
s = Structure(material, grid, solver)
# dump(s.grid)
# display(s.grid.elements[1].cv.N); println()
# display(s.grid.elements[1].cv.dNdξ); println()

# add displacement boundary conditions
node_ids = find_nodes(s, (-1,-1,-1), (1e-3, 1e2, 1e2))
cdofs = [1,2,3]
add_bc!(s, node_ids, cdofs, t -> (0,0,0))

# add external loads
node_ids = find_nodes(s, (2-0.1, 1-0.1, -1e2), (1e2, 1e2, 1e2))
add_force!(s, node_ids, (x,t) -> (0,-1e7,0))

# save data
N = 1000
save(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/beam3d_static/structure_"*string(N))

# solve
solve!(s)

# save final data
# save(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/beam3d_static/structure_"*string(N+1))