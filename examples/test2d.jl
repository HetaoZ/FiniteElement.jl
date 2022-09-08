include("../src/FiniteElement.jl")
using .FiniteElement

nel = (2,1)
grid = RectangularGrid{2,Quadrilateral}((0,0), (2,1), nel)

# dump(grid)

E = 1e9
ν = 0.3
ρ₀ = 1
material = LinearElasticity(2,E,ν,ρ₀)
# display(material);println()

solver = NewtonRaphsonSolver

s = Structure(material, grid, solver)

# 约束边界
node_ids = find_nodes(s, (-1,-1), (1e-6,1e2))
function bc!(node::Node)
    node.d = Vec(0,0)
end

# add_dirichlet!(s, node_ids, bc!)

# 施加载荷
# node_ids = find_nodes(s, (2-0.1,0.1), (2+0.1,1e2))
# add_node_force!(s, node_ids, t->Vec(0,-1e3))

# 导出
# dump(s)
N = 1000000
save_vtk(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/test2d/structure_"*string(N))


# 求解
advance!(s)