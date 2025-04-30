using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 1.0
nf = (4,4)
F1,V1 = hexagonmesh(r,nf; weave=0.0)

## Visualize mesh
markersize = 25
linewidth = 2

fig = Figure(size = (1200,800))  
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "hexagonmesh")
hp1 = poly!(ax1, GeometryBasics.Mesh(V1,F1), color=:white,
strokecolor=:black, strokewidth=3,
transparency=false,shading = FastShading)
hp2 = scatter!(ax1,V1,markersize=markersize,color=:black)
fig