using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

plateDim1 = [20.0,24.0]
pointSpacing1 = 2.0

orientation1 = :up
F1,V1 = triplate(plateDim1,pointSpacing1; orientation=orientation1)

# Visualization
GLMakie.closeall()

fig = Figure(size=(1200,800))
ax1 = AxisGeom(fig[1, 1], title = "Triangulated mesh plate")
hp2 = meshplot!(ax1, F1, V1)
# normalplot(ax1,F1,V1)
fig
