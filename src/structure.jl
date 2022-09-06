

function Structure(material::AbstractMaterial, grid::Grid{dim,T}, solver::AbstractSolver) where {dim,T}
    return Structure{dim}(material, grid, solver,
     new_solution(solver, getndofs(grid)), 
     AbstractConstrain[], Dict(), true)
end

"根据prototype自动生成网格"
function Structure(material::AbstractMaterial, grid_prototype::AbstractGridPrototype, solver::AbstractSolver)
    return Structure(material, generate(grid_prototype), solver)
end

"获取指定空间范围内包含的结点编号向量"
function find_nodes(s::Structure, start, stop)
    return find_nodes(s.grid, start, stop)
end

"""
Fetch data from elements. This differs from assemble_elem_field.
"""
function fetch_data(s::Structure{dim}, field) where dim
    f = getfield(s.grid.nodes[1], field)
    d = Array{eltype(f)}(undef, dim, getnnodes(s.grid))
    for node in s.grid.nodes 
        d[:,node.id] = getfield(node, field)
    end
    return d    
end