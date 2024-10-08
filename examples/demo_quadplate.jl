using Comodo
using GLMakie
using GeometryBasics

plateDim1 = [20.0,24.0]
plateElem1 = [11,16]
orientation1 = :up
F1,V1 = quadplate(plateDim1,plateElem1; orientation=orientation1)

## Visualize mesh
fig = Figure(size=(1200,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Quadrilateral mesh plate")
hp2 = poly!(ax1,GeometryBasics.Mesh(V1,F1), strokewidth=3,color=:white, shading = FastShading)

fig
