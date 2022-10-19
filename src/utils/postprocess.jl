
# 更多类型请参阅 WriteVTK.jl 的代码
const VTK_CELL_TYPE = Dict(
    Line=>VTKCellTypes.VTK_LINE, 
    Triangle=>VTKCellTypes.VTK_TRIANGLE, 
    Quadrilateral=>VTKCellTypes.VTK_QUAD,
    Tetrahedron=>VTKCellTypes.VTK_TETRA,
    Hexahedron=>VTKCellTypes.VTK_HEXAHEDRON
)

function save(s::Structure, filename::String) 
    save(s, ("x0","x","d","u","a","σᵥ","σ"), (:x0,:x,:d,:u,:a,:σᵥ,:σ), filename)
end

function save(s::Structure, datanames::Tuple{Vararg{String}}, fields::Tuple{Vararg{Symbol}}, filename::String) 
    vtkfile = create_vtkfile(s.grid, filename)
    for i in eachindex(datanames)
        if fields[i] == :σ
            vtkfile[datanames[i]] = fetch_σ_data(s)
        elseif fields[i] == :σᵥ
            vtkfile[datanames[i]] = fetch_σᵥ_data(s)
        else
            vtkfile[datanames[i]] = fetch_data(s, fields[i])
        end
    end
    outfiles = vtk_save(vtkfile)
    # println("       saved to ",outfiles[1])
end

function create_vtkfile(grid::Grid, fname)
    points, cells = pts_cells(grid)
    f = vtk_grid(fname, points, cells)
    return f
end

function pts_cells(grid::Grid{dim,T}) where {dim,T}

    pts = Matrix{Float64}(undef, dim, getnnodes(grid))
    cells = Vector{MeshCell}(undef, getnelems(grid))

    for eid in eachindex(grid.elements)
        e = grid.elements[eid]
        
        for k in eachindex(e.connection)
            node = grid.nodes[e.connection[k]]
            pts[:, node.id] = node.x0
        end
        cells[eid] = MeshCell(VTK_CELL_TYPE[T], e.connection)
    end
    return pts, cells
end

"""
Fetch data from nodes.
"""
function fetch_data(s::Structure{dim}, field) where dim
    f = getfield(s.grid.nodes[1], field)
    d = Array{eltype(f)}(undef, dim, getnnodes(s.grid))
    for node in s.grid.nodes 
        d[:,node.id] = getfield(node, field)
    end
    return d    
end

"""
Fetch data from cells (converted to node data).
"""
function fetch_σ_data(s::Structure{dim}) where dim
    d = zeros(Float64, dim, getnnodes(s.grid))
    cell_count = zeros(Int, getnnodes(s.grid))
    # if typeof(s.grid.elements[1]) <: SimplexElement
        for (elem_id, elem) in enumerate(s.grid.elements)
            for i in eachindex(elem.connection)
                node_id = elem.connection[i]
                cell_count[node_id] += 1
                for axis in 1:dim
                    d[axis, node_id] += mean([state.σ[axis,axis] for state in s.states[elem_id]])
                end
            end
        end
        for axis in 1:dim
            d[axis, :] = d[axis,:] ./ cell_count
        end

    # else
    #     error("undef elemtype")
    # end
    return d    
end

"""
Fetch data from cells (converted to node data).
"""
function fetch_σᵥ_data(s::Structure{dim}) where dim
    d = zeros(Float64, getnnodes(s.grid))
    cell_count = zeros(Int, getnnodes(s.grid))
    # if typeof(s.grid.elements[1]) <: SimplexElement
        for (elem_id, elem) in enumerate(s.grid.elements)
            for i in eachindex(elem.connection)
                node_id = elem.connection[i]
                cell_count[node_id] += 1
                d[node_id] += mean([state.σᵥ for state in s.states[elem_id]])
            end
        end
        d = reshape(d ./ cell_count, (1,getnnodes(s.grid)))
    # else
    #     error("undef elemtype")
    # end
    return d    
end


import Base.show

function show(s::Structure{dim}) where {dim} 
    println(dim,"-D Structure")
    println("--  ", length(s.grid.elements) , " elements of type ", eltype(s.grid.elements))
end

