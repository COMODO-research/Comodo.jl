using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 1
F,V = dodecahedron(r)

## Visualize mesh
GLMakie.closeall()

markersize = 25
strokewidth = 2 
strokecolor = :black

fig = Figure(size = (800,800))
ax1 = AxisGeom(fig[1, 1], title = "Dodecahedron")
hp1 = meshplot!(ax1, F, V, strokewidth=strokewidth, strokecolor=strokecolor)
hp2 = scatter!(ax1, V,markersize=markersize, color=:red)
hp3 = normalplot(ax1, F, V; color = :green)
fig