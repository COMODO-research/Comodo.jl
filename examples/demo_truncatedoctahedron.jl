using Comodo
using GLMakie
using GeometryBasics
using Rotations

# Creating faces and vertices for a truncated octahedron
w = 1.0 
F,V = truncatedoctahedron(w)

## Visualize mesh
markersize = 25
strokewidth = 2 
strokecolor = :black
fig = Figure(size = (800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "truncatedoctahedron")

hp1 = poly!(ax1, GeometryBasics.Mesh(V,F[1]), color=:white,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)
hp2 = poly!(ax1, GeometryBasics.Mesh(V,F[2]), color=:white,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)
hp3 = scatter!(ax1, V,markersize=markersize,color=:red)
hp4 = normalplot(ax1,F[1],V; color = :green)
hp5 = normalplot(ax1,F[2],V; color = :green)

fig