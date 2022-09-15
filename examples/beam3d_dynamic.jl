include("../src/FiniteElement.jl")
using .FiniteElement

# define a 2D rectangular grid
nel = (4,2,2)
grid = RectangularGrid{3, Hexahedron}((0,0,0), (2,1,1), nel)


# define a material
E = 1e9
ν = 0.3
ρ₀ = 1
material = LinearElasticity(3,E,ν,ρ₀)

# define a solver for static analysis
solver = NewmarkSolver

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
add_force!(s, node_ids, (x,t) -> (0,-1e6,0))
"需要想办法添加指定坐标的点力、面力、体积力"

# save data
N = 1000
save(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/beam3d_dynamic/structure_"*string(N))

# solve
t = 0
for i in 1:10
    global t
    println("i = ", i, "  t = ", t)

    Δt = time_step!(s)
    solve!(s, Δt, t)
    t += Δt

    # save final data
    save(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/beam3d_dynamic/structure_"*string(N+i))

    
end