using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 1.0
nf = (9,6)
V = hexagongrid(r,nf; weave=0.0)

## Visualize mesh
markersize = 25
linewidth = 2

fig = Figure(size = (1200,800))  
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "hexagonmesh")
hp1 = scatter!(ax1,V,markersize=markersize,color=:black)
fig