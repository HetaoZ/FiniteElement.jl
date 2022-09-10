tuplezeros(n) = Tuple(zeros(n))
tensorzeros(n) = Vec(tuplezeros(n))

function vonMises(σ, dim)
    if dim == 3
        s = dev(σ)
        return sqrt(3.0/2.0 * s ⊡ s)
    else
        error("undef")
    end
end

@inline function add_cartesian(id::CartesianIndex, axis::Int, n::Int)
    return CartesianIndex(ntuple(k -> k == axis ? id[k]+n : id[k], length(id)))
end

function polygon_area(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    if n < 3 
        return 0.
    end
    s = x[n]*(y[1] - y[n-1]) + x[1]*(y[2] - y[n])
    for i = 2:n-1
        s += x[i]*(y[i+1] - y[i-1])
    end
    return s * 0.5
end

"""
注意这里的结点必须按如下顺序排列：

第一层：

P3--P4

|    |

P1--P2

第二层：

P7--P8

|    |

P5--P6

计算体积时要考虑到六面体的每个面都是双线性插值曲面，因此每个面不能简单划分为两个三角形，而应该划分为四个三角形。

尤其要注意的是，计算体积时必须确保所有facet的法向都指向六面体外侧。
"""
function hexahedron_volume(P1, P2, P3, P4, P5, P6, P7, P8)
    # 每个面上的四个点的排列顺序是固定的，使得每个facet的法向都指向六面体外侧
    V1 = four_facets_volume(P1, P5, P7, P3)
    V2 = four_facets_volume(P2, P4, P8, P6)

    V3 = four_facets_volume(P1, P2, P6, P5)
    V4 = four_facets_volume(P4, P3, P7, P8)

    V5 = four_facets_volume(P1, P2, P4, P3)
    V6 = four_facets_volume(P5, P6, P8, P7)

    V = abs(V1+V2+V3+V4+V5+V6)
    return V
end

"如果相对于原点，4个顶点逆时针排列（任意三个位矢组成右手系）则体积为正，否则为负。"
function four_facets_volume(P1, P2, P3, P4)
    C = (P1+P2+P3+P4) / 4
    V = facet_volume(P1,P2,C) + facet_volume(P2,P3,C) + facet_volume(P3,P4,C) + facet_volume(P4,P1,C)
    return V
end

"""
三维空间中一个三角形面片facet与原点构成的四面体体积，带有正负符号。如果相对于原点，三个顶点逆时针排列（三个位矢组成右手系）则体积为正，否则为负。
"""
function facet_volume(P1, P2, P3)
    return mixed_product(P1,P2,P3) / 6
end

function tetrahedron_volume(P1, P2, P3, P4)
    return abs(mixed_product(P2-P1, P3-P1, P4-P1)) / 6
end

"Mixed product of three 3-vectors (a, b, c) = dot(a, cross(b,c)) where dot() and cross() are imported from LinearAlgebra"
function mixed_product(a::Vector{Float64},b::Vector{Float64},c::Vector{Float64})
    return dot(a, cross(b, c))
end

@inline function between(a, lo, hi)
    return all(lo .< a .< hi)
end

@inline function betweeneq(a, lo, hi)
    return all(lo .≤ a .≤ hi)
end