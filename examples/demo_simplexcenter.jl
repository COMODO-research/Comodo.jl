using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 10.0 # radius
n = 1 # Number of refinement iterations
F1,V1 = geosphere(n,r)

plateDim = (24,24)
plateElem = (5,5)

F2,V2 = quadplate(plateDim,plateElem)

V1C = simplexcenter(F1,V1)
V2C = simplexcenter(F2,V2)

## Visualization parameters
stroke_width = 2
marker_size = 15 

## Visualize mesh
GLMakie.closeall()

fig = Figure(size=(800,800))

ax1 = AxisGeom(fig[1, 1], title = "Triangle centre points")
hp1 = meshplot!(ax1, F1, V1, strokewidth=stroke_width)
hs2 = scatter!(ax1, V1C, markersize=marker_size, color=:blue)

ax2 = AxisGeom(fig[1, 2], title = "Quad centre points")
hp2 = meshplot!(ax2, F2, V2, strokewidth=stroke_width)
hs2 = scatter!(ax2, V2C, markersize=marker_size, color=:blue)

fig
