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

@inline function between(a, lo, hi)
    return all(lo .< a .< hi)
end

@inline function betweeneq(a, lo, hi)
    return all(lo .≤ a .≤ hi)
end