include("../src/FiniteElement.jl")
using .FiniteElement

# define a 2D rectangular grid
grid = read_grid("/home/hetao/Downloads/models-master-Waveguides/models-master-Waveguides/Waveguides/waveguide3d.msh")


# define a material
E = 1e3
ν = 0.3
ρ₀ = 1
H = E/20
σ₀ = 1
# material = LinearElasticity(3,E,ν,ρ₀)
material = J2Plasticity(3,E,ν,σ₀,H,ρ₀)

# define a solver for static analysis
# solver = ExplicitSolver
solver = NewmarkSolver

# create a structure
s = Structure(material, grid, solver)
show(s)


λ = 1e9 # frequency
T = 1/λ

# # add displacement boundary conditions
node_ids = find_nodes(s, (0.254358,0.482822,0.025), (0.290958,0.461692,0.025), (0.254358,0.482822,-0.025))
node_ids = append!(node_ids, find_nodes(s, (0.290958,-0.461692,0.025), (0.254358,-0.482822,0.025), (0.290958,-0.461692,-0.025)))
cdofs = [1,2,3]
add_disp!(s, node_ids, cdofs, (x,t) -> (0,0,0))

node_ids = find_nodes(s, (-0.545315,-0.0211309,0.025), (-0.545315,-0.0211309,-0.025), (-0.545315,0.0211309,-0.025))
cdofs = [1,2,3]

add_disp!(s, node_ids, cdofs, (x,t) -> (5e-2*sin(t/T*2π),0,0))

# # add external loads

# add_force!(s, node_ids, (x,t) ->  (1e11*sin(t/T), 0, 0) )
# add_force!(s, node_ids, (x,t) -> t < 7*T ? (1e11*sin(t/T), 0, 0) : (0,0,0))
# add_force!(s, node_ids, (x,t) -> t < 5*T ? (5e10, 0, 0) : (0,0,0))
# "需要想办法添加指定坐标的点力、面力、体积力"

# display(s.constrains);println()

# save data
N = 1000000
save(s, "../../out/waveguide3d_dynamic/structure_"*string(N))

# # solve
t = 0
for i in 1:200
    global t
    global λ
    global T

    Δt = T * 0.01
    # Δt = time_step(s)
    solve!(s, Δt, t)
    t += Δt

    if i%1 == 0
        println("i = ", i, "  t = ", t, "  Δt = ",Δt,"  d = ", 1e-3*sin(t/T*2π))
        save(s, "../../out/waveguide3d_dynamic/structure_"*string(N+i))
    end
    
end