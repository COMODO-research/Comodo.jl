using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

# Creating faces and vertices for a truncated octahedron
w = 1.0 
F,V = truncatedoctahedron(w)

## Visualize mesh
GLMakie.closeall()

markersize = 25
strokewidth = 2 
strokecolor = :black
fig = Figure(size = (800,800))

ax1 = AxisGeom(fig[1, 1], title = "truncatedoctahedron")

hp1 = meshplot!(ax1, F[1], V, strokewidth=strokewidth, strokecolor=strokecolor)
hp2 = meshplot!(ax1, F[2], V, strokewidth=strokewidth, strokecolor=strokecolor)
hp3 = scatter!(ax1, V, markersize=markersize, color=:red)
hp4 = normalplot(ax1, F[1], V; color = :green)
hp5 = normalplot(ax1, F[2], V; color = :green)

fig