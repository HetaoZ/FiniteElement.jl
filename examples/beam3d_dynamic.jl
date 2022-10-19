include("../src/FiniteElement.jl")
using .FiniteElement

# define a 2D rectangular grid
nel = (10,1,1) .* 1
# grid = RectangularGrid{3, Hexahedron}((0,0,0), (2,1,1), nel)
grid = RectangularGrid{3, Tetrahedron}((0,0,0), (10,1,1), nel)


# define a material
E = 220e9
ν = 0.3
ρ₀ = 7.6e3
H = E/20
σ₀ = 200e6
material = LinearElasticity(3,E,ν,ρ₀)
# material = J2Plasticity(3,E,ν,σ₀,H,ρ₀)

# define a solver for static analysis
# solver = ExplicitSolver
solver = NewmarkSolver

# create a structure
s = Structure(material, grid, solver)
# dump(s.grid)
# display(s.grid.elements[1].cv.N); println()
# display(s.grid.elements[1].cv.dNdξ); println()

# add displacement boundary conditions
node_ids = find_nodes(s, (-1,-1,-1), (1e-3, 1e2, 1e2))
node_ids = append!(node_ids, find_nodes(s, (10-1e-3, 0-1e-3, 0-1e-3), (10+1e-3, 1+1e-3, 1+1e-3)))
cdofs = [1,2,3]
add_velocity!(s, node_ids, cdofs, t -> (0,0,0))

# add external loads
node_ids = find_nodes(s, (5-1e-3, 0.0-1e-3, 1-1e-3), (5+1e-3, 1.0+1e-3, 1+1e-3))
# add_force!(s, node_ids, (x,t) -> t < 100e-6 ? (1e7,0,0) : (0,0,0))
add_force!(s, node_ids, (x,t) -> (0,0,-1e9))
"需要想办法添加指定坐标的点力、面力、体积力"

# display(s.constrains);println()
# display(s.grid.nodes);println()

# save data
N = 1000000
save(s, "../../out/beam3d_dynamic/structure_"*string(N))

# solve
t = 0
for i in 1:1000
    global t

    Δt = time_step(s)
    solve!(s, Δt, t)
    t += Δt

    if i%20 == 0
        println("i = ", i, "  t = ", t)
        save(s, "../../out/beam3d_dynamic/structure_"*string(N+i))
    end
    
end