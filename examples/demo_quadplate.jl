using Comodo
using GLMakie
using GeometryBasics

plateDim = [20.0,24.0]
plateElem = [11,16]

F,V = quadplate(plateDim,plateElem)
N = facenormal(F,V)


## Visualize mesh
fig = Figure(size=(1000,1000))

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A quadrilateral mesh of a plate",azimuth=-pi/2,elevation=pi/2)
hp2=poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white, shading = FastShading)
# hp2 = scatter!(ax1, V,markersize=15,color=:orange)
# normalplot(ax1,F,V)
fig
