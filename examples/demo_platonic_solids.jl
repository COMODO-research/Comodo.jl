using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 1.0 #radius
F1,V1 = platonicsolid(1,r)
F2,V2 = platonicsolid(2,r)
F3,V3 = platonicsolid(3,r)
F4,V4 = platonicsolid(4,r)
F5,V5 = platonicsolid(5,r)

## Visualize mesh
GLMakie.closeall()

fig = Figure(size = (1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Tetrahedron")
hp1 = meshplot!(ax1, F1, V1, strokewidth=3,color=:red)

ax2 = AxisGeom(fig[1, 2], title = "Hexahedron (cube)")
hp2 = meshplot!(ax2, F2, V2, strokewidth=3,color=:green)

ax3 = AxisGeom(fig[1, 3], title = "Octahedron")
hp3 = meshplot!(ax3, F3, V3, strokewidth=3,color=:blue)

ax4 = AxisGeom(fig[2, 1], title = "Icosahedron")
hp4 = meshplot!(ax4, F4, V4, strokewidth=3,color=:yellow)

ax5 = AxisGeom(fig[2, 2], title = "Dodecahedron")
hp5 = meshplot!(ax5, F5, V5, strokewidth=3,color=:magenta)

fig