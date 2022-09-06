############
getlowerdim(::Lagrange{dim,shape,order}) where {dim,shape,order} = Lagrange{dim-1,shape,order}()
getlowerorder(::Lagrange{dim,shape,order}) where {dim,shape,order} = Lagrange{dim,shape,order-1}()
getlowerorder(::Lagrange{dim,shape,1}) where {dim,shape} = DiscontinuousLagrange{dim,shape,0}()

##################################
# Lagrange dim 1 RefCube order 1 #
##################################
getnbasefunctions(::Lagrange{1,RefCube,1}) = 2
nvertexdofs(::Lagrange{1,RefCube,1}) = 1

faces(::Lagrange{1,RefCube,1}) = ((1,), (2,))

function reference_coordinates(::Lagrange{1,RefCube,1})
    return [Vec{1, Float64}((-1.0,)),
            Vec{1, Float64}(( 1.0,))]
end

function value(ip::Lagrange{1,RefCube,1}, i::Int, ξ::Vec{1})
    ξ_x = ξ[1]
    i == 1 && return (1 - ξ_x) * 0.5
    i == 2 && return (1 + ξ_x) * 0.5
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

##################################
# Lagrange dim 1 RefCube order 2 #
##################################
getnbasefunctions(::Lagrange{1,RefCube,2}) = 3
nvertexdofs(::Lagrange{1,RefCube,2}) = 1
ncelldofs(::Lagrange{1,RefCube,2}) = 1

faces(::Lagrange{1,RefCube,2}) = ((1,), (2,))

function reference_coordinates(::Lagrange{1,RefCube,2})
    return [Vec{1, Float64}((-1.0,)),
            Vec{1, Float64}(( 1.0,)),
            Vec{1, Float64}(( 0.0,))]
end

function value(ip::Lagrange{1,RefCube,2}, i::Int, ξ::Vec{1})
    ξ_x = ξ[1]
    i == 1 && return ξ_x * (ξ_x - 1) * 0.5
    i == 2 && return ξ_x * (ξ_x + 1) * 0.5
    i == 3 && return 1 - ξ_x^2
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

##################################
# Lagrange dim 2 RefCube order 1 #
##################################
getnbasefunctions(::Lagrange{2,RefCube,1}) = 4
nvertexdofs(::Lagrange{2,RefCube,1}) = 1

faces(::Lagrange{2,RefCube,1}) = ((1,2), (2,3), (3,4), (4,1))

function reference_coordinates(::Lagrange{2,RefCube,1})
    return [Vec{2, Float64}((-1.0, -1.0)),
            Vec{2, Float64}(( 1.0, -1.0)),
            Vec{2, Float64}(( 1.0,  1.0,)),
            Vec{2, Float64}((-1.0,  1.0,))]
end

function value(ip::Lagrange{2,RefCube,1}, i::Int, ξ::Vec{2})
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    i == 1 && return (1 - ξ_x) * (1 - ξ_y) * 0.25
    i == 2 && return (1 + ξ_x) * (1 - ξ_y) * 0.25
    i == 3 && return (1 + ξ_x) * (1 + ξ_y) * 0.25
    i == 4 && return (1 - ξ_x) * (1 + ξ_y) * 0.25
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

##################################
# Lagrange dim 2 RefCube order 2 #
##################################
getnbasefunctions(::Lagrange{2,RefCube,2}) = 9
nvertexdofs(::Lagrange{2,RefCube,2}) = 1
nfacedofs(::Lagrange{2,RefCube,2}) = 1
ncelldofs(::Lagrange{2,RefCube,2}) = 1

faces(::Lagrange{2,RefCube,2}) = ((1,2,5), (2,3,6), (3,4,7), (4,1,8))

function reference_coordinates(::Lagrange{2,RefCube,2})
    return [Vec{2, Float64}((-1.0, -1.0)),
            Vec{2, Float64}(( 1.0, -1.0)),
            Vec{2, Float64}(( 1.0,  1.0)),
            Vec{2, Float64}((-1.0,  1.0)),
            Vec{2, Float64}(( 0.0, -1.0)),
            Vec{2, Float64}(( 1.0,  0.0)),
            Vec{2, Float64}(( 0.0,  1.0)),
            Vec{2, Float64}((-1.0,  0.0)),
            Vec{2, Float64}(( 0.0,  0.0))]
end

function value(ip::Lagrange{2,RefCube,2}, i::Int, ξ::Vec{2})
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    i == 1 && return (ξ_x^2 - ξ_x) * (ξ_y^2 - ξ_y) * 0.25
    i == 2 && return (ξ_x^2 + ξ_x) * (ξ_y^2 - ξ_y) * 0.25
    i == 3 && return (ξ_x^2 + ξ_x) * (ξ_y^2 + ξ_y) * 0.25
    i == 4 && return (ξ_x^2 - ξ_x) * (ξ_y^2 + ξ_y) * 0.25
    i == 5 && return (1 - ξ_x^2) * (ξ_y^2 - ξ_y) * 0.5
    i == 6 && return (ξ_x^2 + ξ_x) * (1 - ξ_y^2) * 0.5
    i == 7 && return (1 - ξ_x^2) * (ξ_y^2 + ξ_y) * 0.5
    i == 8 && return (ξ_x^2 - ξ_x) * (1 - ξ_y^2) * 0.5
    i == 9 && return (1 - ξ_x^2) * (1 - ξ_y^2)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

#########################################
# Lagrange dim 2 RefTetrahedron order 1 #
#########################################
getnbasefunctions(::Lagrange{2,RefTetrahedron,1}) = 3
getlowerdim(::Lagrange{2, RefTetrahedron, order}) where {order} = Lagrange{1, RefCube, order}()
nvertexdofs(::Lagrange{2,RefTetrahedron,1}) = 1

vertices(::Lagrange{2,RefTetrahedron,1}) = (1,2,3)
faces(::Lagrange{2,RefTetrahedron,1}) = ((1,2), (2,3), (3,1))

function reference_coordinates(::Lagrange{2,RefTetrahedron,1})
    return [Vec{2, Float64}((1.0, 0.0)),
            Vec{2, Float64}((0.0, 1.0)),
            Vec{2, Float64}((0.0, 0.0))]
end

function value(ip::Lagrange{2,RefTetrahedron,1}, i::Int, ξ::Vec{2})
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    i == 1 && return ξ_x
    i == 2 && return ξ_y
    i == 3 && return 1. - ξ_x - ξ_y
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

#########################################
# Lagrange dim 2 RefTetrahedron order 2 #
#########################################
getnbasefunctions(::Lagrange{2,RefTetrahedron,2}) = 6
nvertexdofs(::Lagrange{2,RefTetrahedron,2}) = 1
nfacedofs(::Lagrange{2,RefTetrahedron,2}) = 1

vertices(::Lagrange{2,RefTetrahedron,2}) = (1,2,3)
faces(::Lagrange{2,RefTetrahedron,2}) = ((1,2,4), (2,3,5), (3,1,6))

function reference_coordinates(::Lagrange{2,RefTetrahedron,2})
    return [Vec{2, Float64}((1.0, 0.0)),
            Vec{2, Float64}((0.0, 1.0)),
            Vec{2, Float64}((0.0, 0.0)),
            Vec{2, Float64}((0.5, 0.5)),
            Vec{2, Float64}((0.0, 0.5)),
            Vec{2, Float64}((0.5, 0.0))]
end

function value(ip::Lagrange{2,RefTetrahedron,2}, i::Int, ξ::Vec{2})
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    γ = 1. - ξ_x - ξ_y
    i == 1 && return ξ_x * (2ξ_x - 1)
    i == 2 && return ξ_y * (2ξ_y - 1)
    i == 3 && return γ * (2γ - 1)
    i == 4 && return 4ξ_x * ξ_y
    i == 5 && return 4ξ_y * γ
    i == 6 && return 4ξ_x * γ
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

#########################################
# Lagrange dim 3 RefTetrahedron order 1 #
#########################################
getnbasefunctions(::Lagrange{3,RefTetrahedron,1}) = 4
nvertexdofs(::Lagrange{3,RefTetrahedron,1}) = 1

faces(::Lagrange{3,RefTetrahedron,1}) = ((1,2,3), (1,2,4), (2,3,4), (1,4,3))
edges(::Lagrange{3,RefTetrahedron,1}) = ((1,2), (2,3), (3,1), (1,4), (2,4), (3,4))

function reference_coordinates(::Lagrange{3,RefTetrahedron,1})
    return [Vec{3, Float64}((0.0, 0.0, 0.0)),
            Vec{3, Float64}((1.0, 0.0, 0.0)),
            Vec{3, Float64}((0.0, 1.0, 0.0)),
            Vec{3, Float64}((0.0, 0.0, 1.0))]
end

function value(ip::Lagrange{3,RefTetrahedron,1}, i::Int, ξ::Vec{3})
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    ξ_z = ξ[3]
    i == 1 && return 1.0 - ξ_x - ξ_y - ξ_z
    i == 2 && return ξ_x
    i == 3 && return ξ_y
    i == 4 && return ξ_z
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

#########################################
# Lagrange dim 3 RefTetrahedron order 2 #
#########################################
getnbasefunctions(::Lagrange{3,RefTetrahedron,2}) = 10
nvertexdofs(::Lagrange{3,RefTetrahedron,2}) = 1
nedgedofs(::Lagrange{3,RefTetrahedron,2}) = 1

faces(::Lagrange{3,RefTetrahedron,2}) = ((1,2,3,5,6,7), (1,2,4,5,9,8), (2,3,4,6,10,9), (1,4,3,8,10,7))
edges(::Lagrange{3,RefTetrahedron,2}) = ((1,5,2), (2,6,3), (3,7,1), (1,8,4), (2,9,4), (3,10,4))

function reference_coordinates(::Lagrange{3,RefTetrahedron,2})
    return [Vec{3, Float64}((0.0, 0.0, 0.0)),
            Vec{3, Float64}((1.0, 0.0, 0.0)),
            Vec{3, Float64}((0.0, 1.0, 0.0)),
            Vec{3, Float64}((0.0, 0.0, 1.0)),
            Vec{3, Float64}((0.5, 0.0, 0.0)),
            Vec{3, Float64}((0.5, 0.5, 0.0)),
            Vec{3, Float64}((0.0, 0.5, 0.0)),
            Vec{3, Float64}((0.0, 0.0, 0.5)),
            Vec{3, Float64}((0.5, 0.0, 0.5)),
            Vec{3, Float64}((0.0, 0.5, 0.5))]
end

# http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
# http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch10.d/AFEM.Ch10.pdf
function value(ip::Lagrange{3,RefTetrahedron,2}, i::Int, ξ::Vec{3})
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    ξ_z = ξ[3]
    i == 1  && return (-2 * ξ_x - 2 * ξ_y - 2 * ξ_z + 1) * (-ξ_x - ξ_y - ξ_z + 1)
    i == 2  && return ξ_x * (2 * ξ_x - 1)
    i == 3  && return ξ_y * (2 * ξ_y - 1)
    i == 4  && return ξ_z * (2 * ξ_z - 1)
    i == 5  && return ξ_x * (-4 * ξ_x - 4 * ξ_y - 4 * ξ_z + 4)
    i == 6  && return 4 * ξ_x * ξ_y
    i == 7  && return 4 * ξ_y * (-ξ_x - ξ_y - ξ_z + 1)
    i == 8  && return ξ_z * (-4 * ξ_x - 4 * ξ_y - 4 * ξ_z + 4)
    i == 9  && return 4 * ξ_x * ξ_z
    i == 10 && return 4 * ξ_y * ξ_z
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end

##################################
# Lagrange dim 3 RefCube order 1 #
##################################
getnbasefunctions(::Lagrange{3,RefCube,1}) = 8
nvertexdofs(::Lagrange{3,RefCube,1}) = 1

faces(::Lagrange{3,RefCube,1}) = ((1,4,3,2), (1,2,6,5), (2,3,7,6), (3,4,8,7), (1,5,8,4), (5,6,7,8))
edges(::Lagrange{3,RefCube,1}) = ((1,2), (2,3), (3,4), (4,1), (1,5), (2,6), (3,7), (4,8), (5,6), (6,7), (7,8), (8,5))

function reference_coordinates(::Lagrange{3,RefCube,1})
    return [Vec{3, Float64}((-1.0, -1.0, -1.0)),
            Vec{3, Float64}(( 1.0, -1.0, -1.0)),
            Vec{3, Float64}(( 1.0,  1.0, -1.0)),
            Vec{3, Float64}((-1.0,  1.0, -1.0)),
            Vec{3, Float64}((-1.0, -1.0,  1.0)),
            Vec{3, Float64}(( 1.0, -1.0,  1.0)),
            Vec{3, Float64}(( 1.0,  1.0,  1.0)),
            Vec{3, Float64}((-1.0,  1.0,  1.0))]
end

function value(ip::Lagrange{3,RefCube,1}, i::Int, ξ::Vec{3})
    ξ_x = ξ[1]
    ξ_y = ξ[2]
    ξ_z = ξ[3]
    i == 1 && return 0.125(1 - ξ_x) * (1 - ξ_y) * (1 - ξ_z)
    i == 2 && return 0.125(1 + ξ_x) * (1 - ξ_y) * (1 - ξ_z)
    i == 3 && return 0.125(1 + ξ_x) * (1 + ξ_y) * (1 - ξ_z)
    i == 4 && return 0.125(1 - ξ_x) * (1 + ξ_y) * (1 - ξ_z)
    i == 5 && return 0.125(1 - ξ_x) * (1 - ξ_y) * (1 + ξ_z)
    i == 6 && return 0.125(1 + ξ_x) * (1 - ξ_y) * (1 + ξ_z)
    i == 7 && return 0.125(1 + ξ_x) * (1 + ξ_y) * (1 + ξ_z)
    i == 8 && return 0.125(1 - ξ_x) * (1 + ξ_y) * (1 + ξ_z)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end


##################################
# Lagrange dim 3 RefCube order 2 #
##################################
# Based on vtkTriQuadraticHexahedron (see https://kitware.github.io/vtk-examples/site/Cxx/GeometricObjects/IsoparametricCellsDemo/)
getnbasefunctions(::Lagrange{3,RefCube,2}) = 27
nvertexdofs(::Lagrange{3,RefCube,2}) = 1
nedgedofs(::Lagrange{3,RefCube,2}) = 1
nfacedofs(::Lagrange{3,RefCube,2}) = 1
ncelldofs(::Lagrange{3,RefCube,2}) = 1

faces(::Lagrange{3,RefCube,2}) = ((1,2,6,5, 9,18,13,17, 23), (2,3,7,6, 10,19,14,18, 22), (3,4,8,7, 11,20,15,19, 24), (1,5,8,4, 12,17,16,20, 21), (1,4,3,2, 9,10,11,12, 25), (5,6,7,8, 13,14,15,16, 26))
edges(::Lagrange{3,RefCube,2}) = ((1,2, 9), (2,3, 10), (3,4, 11), (4,1, 12), (1,5, 16), (2,6, 19), (3,7, 18), (4,8, 19), (5,6, 13), (6,7, 14), (7,8, 15), (8,5, 16))

function reference_coordinates(::Lagrange{3,RefCube,2})
           # vertex
    return [Vec{3, Float64}((-1.0, -1.0, -1.0)), #  0
            Vec{3, Float64}(( 1.0, -1.0, -1.0)), #  1
            Vec{3, Float64}(( 1.0,  1.0, -1.0)), #  2
            Vec{3, Float64}((-1.0,  1.0, -1.0)), #  3
            Vec{3, Float64}((-1.0, -1.0,  1.0)), #  4
            Vec{3, Float64}(( 1.0, -1.0,  1.0)), #  5
            Vec{3, Float64}(( 1.0,  1.0,  1.0)), #  6
            Vec{3, Float64}((-1.0,  1.0,  1.0)), #  7
            # edge
            Vec{3, Float64}(( 0.0, -1.0, -1.0)), #  8
            Vec{3, Float64}(( 1.0,  0.0, -1.0)),
            Vec{3, Float64}(( 0.0,  1.0, -1.0)),
            Vec{3, Float64}((-1.0,  0.0, -1.0)),
            Vec{3, Float64}(( 0.0, -1.0,  1.0)),
            Vec{3, Float64}(( 1.0,  0.0,  1.0)),
            Vec{3, Float64}(( 0.0,  1.0,  1.0)),
            Vec{3, Float64}((-1.0,  0.0,  1.0)),
            Vec{3, Float64}((-1.0, -1.0,  0.0)),
            Vec{3, Float64}(( 1.0, -1.0,  0.0)),
            Vec{3, Float64}(( 1.0,  1.0,  0.0)),
            Vec{3, Float64}((-1.0,  1.0,  0.0)), # 19
            # face
            Vec{3, Float64}(( 0.0, -1.0,  0.0)), # 20
            Vec{3, Float64}(( 1.0,  0.0,  0.0)),
            Vec{3, Float64}(( 0.0,  1.0,  0.0)),
            Vec{3, Float64}((-1.0,  0.0,  0.0)),
            Vec{3, Float64}(( 0.0,  0.0, -1.0)),
            Vec{3, Float64}(( 0.0,  0.0,  1.0)), # 25
            # interior
            Vec{3, Float64}((0.0, 0.0, 0.0)),    # 26
            ]
end

function value(ip::Lagrange{3,RefCube,2}, i::Int, ξ::Vec{3, T}) where {T}
    # Some local helpers.
    @inline φ₁(x::T) = -0.5*x*(1-x)
    @inline φ₂(x::T) = (1+x)*(1-x)
    @inline φ₃(x::T) = 0.5*x*(1+x)
    (ξ_x, ξ_y, ξ_z) = ξ
    # vertices
    i == 1 && return φ₁(ξ_x) * φ₁(ξ_y) * φ₁(ξ_z)
    i == 2 && return φ₃(ξ_x) * φ₁(ξ_y) * φ₁(ξ_z)
    i == 3 && return φ₃(ξ_x) * φ₃(ξ_y) * φ₁(ξ_z)
    i == 4 && return φ₁(ξ_x) * φ₃(ξ_y) * φ₁(ξ_z)
    i == 5 && return φ₁(ξ_x) * φ₁(ξ_y) * φ₃(ξ_z)
    i == 6 && return φ₃(ξ_x) * φ₁(ξ_y) * φ₃(ξ_z)
    i == 7 && return φ₃(ξ_x) * φ₃(ξ_y) * φ₃(ξ_z)
    i == 8 && return φ₁(ξ_x) * φ₃(ξ_y) * φ₃(ξ_z)
    # edges
    i ==  9 && return φ₂(ξ_x) * φ₁(ξ_y) * φ₁(ξ_z)
    i == 10 && return φ₃(ξ_x) * φ₂(ξ_y) * φ₁(ξ_z)
    i == 11 && return φ₂(ξ_x) * φ₃(ξ_y) * φ₁(ξ_z)
    i == 12 && return φ₁(ξ_x) * φ₂(ξ_y) * φ₁(ξ_z)
    i == 13 && return φ₂(ξ_x) * φ₁(ξ_y) * φ₃(ξ_z)
    i == 14 && return φ₃(ξ_x) * φ₂(ξ_y) * φ₃(ξ_z)
    i == 15 && return φ₂(ξ_x) * φ₃(ξ_y) * φ₃(ξ_z)
    i == 16 && return φ₁(ξ_x) * φ₂(ξ_y) * φ₃(ξ_z)
    i == 17 && return φ₁(ξ_x) * φ₁(ξ_y) * φ₂(ξ_z)
    i == 18 && return φ₃(ξ_x) * φ₁(ξ_y) * φ₂(ξ_z)
    i == 19 && return φ₃(ξ_x) * φ₃(ξ_y) * φ₂(ξ_z)
    i == 20 && return φ₁(ξ_x) * φ₃(ξ_y) * φ₂(ξ_z)
    # faces
    i == 21 && return φ₂(ξ_x) * φ₁(ξ_y) * φ₂(ξ_z)
    i == 22 && return φ₃(ξ_x) * φ₂(ξ_y) * φ₂(ξ_z)
    i == 23 && return φ₂(ξ_x) * φ₃(ξ_y) * φ₂(ξ_z)
    i == 24 && return φ₁(ξ_x) * φ₂(ξ_y) * φ₂(ξ_z)
    i == 25 && return φ₂(ξ_x) * φ₂(ξ_y) * φ₁(ξ_z)
    i == 26 && return φ₂(ξ_x) * φ₂(ξ_y) * φ₃(ξ_z)
    # interior
    i == 27 && return φ₂(ξ_x) * φ₂(ξ_y) * φ₂(ξ_z)
    throw(ArgumentError("no shape function $i for interpolation $ip"))
end