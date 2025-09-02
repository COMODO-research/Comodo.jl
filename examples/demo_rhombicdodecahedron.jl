using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

# Creating faces and vertices for a rhombic dodecahedron
w = 1.0 
F,V = rhombicdodecahedron(w)

## Visualize mesh
GLMakie.closeall()

markersize = 25
strokewidth = 2 
strokecolor = :black

fig = Figure(size = (800,800))
ax1 = AxisGeom(fig[1, 1], title = "Rhombic dodecahedron")
hp1 = meshplot!(ax1, F, V, strokewidth=strokewidth, strokecolor=strokecolor)
hp2 = scatter!(ax1, V,markersize=markersize, color=:red)
hp3 = normalplot(ax1,F,V; color = :green)
fig