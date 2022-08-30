using Gmsh
gmsh.initialize()
gmsh.model.add("ExtrudeBox")

meshsize = 0.1e-3
points = [(x1,x2) for x1 in (0,1e-3), x2 in (0,1e-3)]

tag = 0
for p in points
    global tag
    tag += 1
    gmsh.model.geo.addPoint(p..., meshsize, tag)
end

lbase = 100
gmsh.model.geo.addLine(1,2,lbase+1)
gmsh.model.geo.addLine(2,4,lbase+2)
gmsh.model.geo.addLine(4,3,lbase+3)
gmsh.model.geo.addLine(3,1,lbase+4)

# 注意线段的方向用符号表示
gmsh.model.geo.addCurveLoop(lbase .+ [1,2,3,4], lbase+5)

fbase = 200
gmsh.model.geo.addPlaneSurface([lbase+5],fbase+1)

# # 将geo的结果映射为entities
# gmsh.model.geo.synchronize()

# # 划分网格
# gmsh.model.mesh.generate()
# # gmsh.model.mesh.recombine()
# # gmsh.model.mesh.refine()

gmsh.model.geo.extrude((2,fbase+1),0,0,1e-3,[4],[2e-3])

function getTags(dim)
    entities = gmsh.model.getEntities(dim)
    tags = []
    for e in entities
        push!(tags, e[2])
    end
    return tags
end

gmsh.model.geo.addSurfaceLoop(getTags(2), fbase+7)

vbase = 300
gmsh.model.geo.addVolume([fbase+7], vbase+1)

# 将geo的结果映射为entities
gmsh.model.geo.synchronize()

gmsh.model.mesh.generate()

# 导出
gmsh.write("box3d.vtk")