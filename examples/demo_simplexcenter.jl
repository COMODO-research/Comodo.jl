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
fig = Figure(size=(800,800))

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Triangle centre points")
hp1=poly!(ax1,V1,F1, strokewidth=stroke_width,color=:white, shading = FastShading, transparency=false)
hs2 = scatter!(ax1, V1C,markersize=marker_size,color=:blue)

ax2=Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Quad centre points")
hp2=poly!(ax2,GeometryBasics.Mesh(V2,F2), strokewidth=stroke_width,color=:white, shading = FastShading, transparency=false)
hs2 = scatter!(ax2, V2C,markersize=marker_size,color=:blue)

fig
