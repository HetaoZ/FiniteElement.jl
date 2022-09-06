
# 更多类型请参阅 WriteVTK.jl 的代码
const VTK_CELL_TYPE = Dict(
    Line=>VTKCellTypes.VTK_LINE, 
    Triangle=>VTKCellTypes.VTK_TRIANGLE, 
    Quadrilateral=>VTKCellTypes.VTK_QUAD,
    Tetrahedron=>VTKCellTypes.VTK_TETRA,
    Hexahedron=>VTKCellTypes.VTK_HEXAHEDRON
)

function save_vtk(s::Structure, datanames::Tuple{Vararg{String}}, fields::Tuple{Vararg{Symbol}}, filename::String) 
    vtkfile = create_vtkfile(s.grid, filename)
    for i in eachindex(datanames)
        vtkfile[datanames[i]] = fetch_data(s, fields[i])
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
