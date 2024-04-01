using Comodo
using GLMakie
using GeometryBasics

plateDim = [20.0,20.0]
plateElem = [10,10]

F,V = quadplate(plateDim,plateElem)

## Visualize mesh
fig = Figure(size=(800,800))

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A quadrilateral mesh of a plate")
hp2=poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white, shading = FastShading)

fig
