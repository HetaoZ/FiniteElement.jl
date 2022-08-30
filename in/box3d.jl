using Gmsh
gmsh.initialize()
gmsh.model.add("box3d")

meshsize = 0.5e-3
points = [(x1,x2,x3) for x1 in (0,1e-3), x2 in (0,1e-3), x3 in (0,1e-3)]

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

gmsh.model.geo.addLine(5,6,lbase+5)
gmsh.model.geo.addLine(6,8,lbase+6)
gmsh.model.geo.addLine(8,7,lbase+7)
gmsh.model.geo.addLine(7,5,lbase+8)

gmsh.model.geo.addLine(1,5,lbase+9)
gmsh.model.geo.addLine(2,6,lbase+10)
gmsh.model.geo.addLine(4,8,lbase+11)
gmsh.model.geo.addLine(3,7,lbase+12)

# 注意线段的方向用符号表示
gmsh.model.geo.addCurveLoop(lbase .+ [1,2,3,4], lbase+13)
gmsh.model.geo.addCurveLoop(lbase .+ [5,6,7,8], lbase+14)
gmsh.model.geo.addCurveLoop(lbase .+ [1,10,-5-2*lbase,-9-2*lbase], lbase+15)
gmsh.model.geo.addCurveLoop(lbase .+ [2,11,-6-2*lbase,-10-2*lbase], lbase+16)
gmsh.model.geo.addCurveLoop(lbase .+ [3,12,-7-2*lbase,-11-2*lbase], lbase+17)
gmsh.model.geo.addCurveLoop(lbase .+ [4,9,-8-2*lbase,-12-2*lbase], lbase+18)

fbase = 200
gmsh.model.geo.addPlaneSurface([lbase+13],fbase+1)
gmsh.model.geo.addPlaneSurface([lbase+14],fbase+2)
gmsh.model.geo.addPlaneSurface([lbase+15],fbase+3)
gmsh.model.geo.addPlaneSurface([lbase+16],fbase+4)
gmsh.model.geo.addPlaneSurface([lbase+17],fbase+5)
gmsh.model.geo.addPlaneSurface([lbase+18],fbase+6)

# 将geo的结果映射为entities
gmsh.model.geo.synchronize()

function getTags(dim)
    entities = gmsh.model.getEntities(dim)
    println(entities)
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

# 划分网格
gmsh.model.mesh.generate()

# recombine只适用于2D
# gmsh.model.mesh.recombine()

# 导出
gmsh.write("box3d.msh")