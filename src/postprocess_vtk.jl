
const VTK_CELL_TYPE = Dict(
    "Tri3"=>VTKCellTypes.VTK_TRIANGLE, 
    "Quad4"=>VTKCellTypes.VTK_QUAD,
    "Tetra4"=>VTKCellTypes.VTK_TETRA,
    "Hex8"=>VTKCellTypes.VTK_HEXAHEDRON
    )

function create_vtkfile(s::PlasticStructure, fname::String)
    points, cells = pts_cells(s)
    f = vtk_grid(fname, points, cells)
    return f
end

function pts_cells(s::PlasticStructure)
    pts = Matrix{Float64}(undef, length(s.grid.nodes[1].x), length(s.grid.nodes))
    ncells = length(s.grid.cells)
    cells = Vector{MeshCell}(undef, ncells)
    for eid = 1:ncells
        e = s.grid.cells[eid]
        for k = 1:length(e.nodes)
            node = s.grid.nodes[e.nodes[k]]
            # 这里的x相当于x0，因为Ferrite不更新node.x
            pts[:, e.nodes[k]] = node.x
        end
        cells[eid] = MeshCell(VTK_CELL_TYPE["Tetra4"], [e.nodes[i] for i=1:length(e.nodes)])
    end
    return pts, cells
end

function save_to_vtk(s::PlasticStructure, datanames::T where T<:Tuple, fields::T where T<:Tuple, fname::String)
    vtkfile = create_vtkfile(s, fname)
    for i in eachindex(datanames)
        vtkfile[datanames[i]] = fetch_data(s, fields[i])
    end
    outfiles = vtk_save(vtkfile)
    # println("       saved to ",outfiles[1])
end
