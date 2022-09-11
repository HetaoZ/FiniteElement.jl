# FiniteElement

A developing package for FEM analysis. DO NOT use it in working environment.

## Installation
```
] add https://github.com/HetaoZ/FiniteElement.jl.git
```

## Examples
```
using FiniteElement

nel = (4,2)
grid = RectangularGrid{2,Quadrilateral}((0,0), (2,1), nel)

# dump(grid)

E = 1e9
ν = 0.3
ρ₀ = 1
material = LinearElasticity(2,E,ν,ρ₀)
# display(material);println()

solver = StaticSolver(1.0)

s = Structure(material, grid, solver)

# 约束边界
node_ids = find_nodes(s, (-1,-1), (1e-6,1e2))
cdofs = [1,2]
add_bc!(s, node_ids, cdofs, t -> (0,0))

# 施加载荷
node_ids = find_nodes(s, (2-0.1, 0.5), (2+0.1, 1e2))
add_force!(s, node_ids, (x,t) -> (0,-1e7))

# 导出
# dump(s)

N = 1000
save(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/test2d/structure_"*string(N))


# 求解
solve!(s)
save(s, ("x0","x","d","u","a"), (:x0,:x,:d,:u,:a), "../../out/test2d/structure_"*string(N+1))
```
