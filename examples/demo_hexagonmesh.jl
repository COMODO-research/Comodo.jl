using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 1.0
nf = (4,4)
F1,V1 = hexagonmesh(r,nf; weave=0.0)

## Visualize mesh
GLMakie.closeall()

markersize = 25
linewidth = 2

fig = Figure(size = (1200,800))  
ax1 = AxisGeom(fig[1, 1], title = "hexagonmesh")
hp1 = meshplot!(ax1, F1, V1, strokewidth=3)
hp2 = scatter!(ax1,V1,markersize=markersize,color=:black)
fig