include("gaussquad_tri_table.jl")
include("gaussquad_tet_table.jl")
include("generate_quadrature.jl")

import Base.Cartesian: @nloops, @nref, @ntuple, @nexprs


getweights(quad_rule::QuadratureRule) = quad_rule.weights
getpoints(quad_rule::QuadratureRule) = quad_rule.points

QuadratureRule{dim,shape}(order::Int) where {dim,shape} = QuadratureRule{dim,shape}(:legendre, order)

# Special case for face integration of 1D problems
function (::Type{QuadratureRule{0, RefCube}})(quad_type::Symbol, order::Int)
    w = Float64[1.0]
    p = Vec{0,Float64}[]
    return QuadratureRule{0,RefCube,Float64}(w,p)
end

# Generate Gauss quadrature rules on cubes by doing an outer product
# over all dimensions
for dim in (1,2,3)
    @eval begin
        function (::Type{QuadratureRule{$dim,RefCube}})(quad_type::Symbol, order::Int)
            if quad_type == :legendre
                p, w = GaussQuadrature.legendre(Float64, order)
            elseif quad_type == :lobatto
                p, w = GaussQuadrature.legendre(Float64, order, GaussQuadrature.both)
            else
                throw(ArgumentError("unsupported quadrature rule"))
            end
            weights = Vector{Float64}(undef, order^($dim))
            points = Vector{Vec{$dim,Float64}}(undef, order^($dim))
            count = 1
            @nloops $dim i j->(1:order) begin
                t = @ntuple $dim q-> p[$(Symbol("i"*"_q"))]
                points[count] = Vec{$dim,Float64}(t)
                weight = 1.0
                @nexprs $dim j->(weight *= w[i_{j}])
                weights[count] = weight
                count += 1
            end
            return QuadratureRule{$dim,RefCube,Float64}(weights, points)
        end
    end
end

for dim in (2, 3)
    @eval begin
        function (::Type{QuadratureRule{$dim, RefTetrahedron}})(quad_type::Symbol, order::Int)
            if $dim == 2 && quad_type == :legendre
                data = _get_gauss_tridata(order)
            elseif $dim == 3 && quad_type == :legendre
                data = _get_gauss_tetdata(order)
            else
                throw(ArgumentError("unsupported quadrature rule"))
            end
            n_points = size(data,1)
            points = Vector{Vec{$dim,Float64}}(undef, n_points)

            for p in axes(data, 1)
                points[p] = Vec{$dim,Float64}(@ntuple $dim i -> data[p, i])
            end
            weights = data[:, $dim + 1]
            QuadratureRule{$dim,RefTetrahedron,Float64}(weights, points)
        end
    end
end

# Special version for face integration of triangles
function (::Type{QuadratureRule{1,RefTetrahedron}})(quad_type::Symbol, order::Int)
    if quad_type == :legendre
        p, weights = GaussQuadrature.legendre(Float64,order)
    elseif quad_type == :lobatto
        p, weights = GaussQuadrature.legendre(Float64, order, GaussQuadrature.both)
    else
        throw(ArgumentError("unsupported quadrature rule"))
    end
    points = Vector{Vec{1,Float64}}(undef, order)
    # Shift interval from (-1,1) to (0,1)
    weights *= 0.5
    p .+= 1.0; p /= 2.0

    for i in 1:length(weights)
        points[i] = Vec{1,Float64}((p[i],))
    end
    return QuadratureRule{1,RefTetrahedron,Float64}(weights, points)
end
