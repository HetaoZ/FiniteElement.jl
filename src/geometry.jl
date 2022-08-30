function volume(cell::Tetrahedron, nodes::Vector{Node{3,T}} where T <: Real, d::Vector{Float64})
    P = ones(Float64, 4, 4)
    for j = 1:4
        noded = getd(d, cell.nodes[j], 3)
        x = nodes[cell.nodes[j]].x
        for i = 1:3
            P[i+1,j] = x[i] + noded[i]
        end
    end
    return abs(det(P)) / 6
end

function volume(cell::Quadrilateral, nodes::Vector{Node{2, T}} where T <: Real, d::Vector{Float64})
    x1 = zeros(4)
    x2 = zeros(4)
    for j = 1:4
        noded = getd(d, cell.nodes[j], 2)
        x = nodes[cell.nodes[j]].x 
        x1[j] = x[1]
        x2[j] = x[2]
    end
    return abs(polygon_area(x1, x2))
end

function polygon_area(x, y)
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

@inline function getd(d, node_id, dim)
    return d[(node_id-1)*dim+1:node_id*dim]
end