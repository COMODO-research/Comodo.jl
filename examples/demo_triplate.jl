using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

plateDim1 = [20.0,24.0]
pointSpacing1 = 2.0

orientation1 = :up
F1,V1 = triplate(plateDim1,pointSpacing1; orientation=orientation1)

## Visualization
linewidth = 1
fig = Figure(size=(1200,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Triangulated mesh plate")
hp2 = poly!(ax1,GeometryBasics.Mesh(V1,F1), strokewidth=linewidth,color=:white, shading = FastShading)
# normalplot(ax1,F1,V1)
fig
