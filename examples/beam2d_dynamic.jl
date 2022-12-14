include("../src/FiniteElement.jl")
using .FiniteElement

# define a 2D rectangular grid
nel = (1,1)
grid = RectangularGrid{2,Triangle}((0,0), (1,1), nel)

# define a material
E = 1e9
ν = 0.3
ρ₀ = 100
material = LinearElasticity(2,E,ν,ρ₀)

# define a solver for dynamic analysis
solver = NewmarkSolver

# create a structure
s = Structure(material, grid, solver)

# add displacement boundary conditions
node_ids = find_nodes(s, (-1,-1), (1e-6,1e2))
cdofs = [1,2]
add_disp!(s, node_ids, cdofs, t -> (0,0))

# add external loads
node_ids = find_nodes(s, (2-0.1, 0.5), (2+0.1, 1e2))
add_force!(s, node_ids, (x,t) -> (0,-1e6))


# save data
N = 1000
save(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/beam2d_dynamic/structure_"*string(N))

# # solve
t = 0
for i in 1:10
    global t
    println("i = ", i, "  t = ", t)

    Δt = time_step(s)
    solve!(s, Δt, t)
    t += Δt

    # save final data
    save(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/beam2d_dynamic/structure_"*string(N+i))

    
end
