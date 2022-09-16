include("../src/FiniteElement.jl")
using .FiniteElement

# define a 2D rectangular grid
nel = (1,1,1)
grid = RectangularGrid{3, Hexahedron}((0,0,0), (1,1,1), nel)


# define a material
E = 1e9
ν = 0.3
ρ₀ = 1
material = LinearElasticity(3,E,ν,ρ₀)

# define a solver for static analysis
solver = StaticSolver(1e-3)

# create a structure
s = Structure(material, grid, solver)

surface = FiniteElement.getsurface!(s)
println(surface.start)
println(surface.stop)

points = (
    (0.5,0.5,0.5),
    (0.1,0.1,0.1),
    (0.2,0.1,0.1),
    (0.3,0.1,0.1),
    (0.4,0.1,0.1),
    (0.5,0.1,0.1),
    (0.6,0.1,0.1),
    (0.7,0.1,0.1),
    (0.9,0.9,0.9),
    (0.5,0.5,1.0),
    (-0.5,0.5,2.0),
    (3.5,0.5,1.0),
    (1.0, 0.5, 0.5)
    )
for p in points
    println(p, " => " , FiniteElement.pinpoly(p, surface))
end