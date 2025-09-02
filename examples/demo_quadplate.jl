using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

plateDim1 = [20.0,24.0]
plateElem1 = [11,16]
orientation1 = :up
F1,V1 = quadplate(plateDim1,plateElem1; orientation=orientation1)

## Visualize mesh
GLMakie.closeall()

fig = Figure(size = (1200,800))
ax1 = AxisGeom(fig[1, 1], title = "Quadrilateral mesh plate")
hp2 = meshplot!(ax1, F1, V1)
fig