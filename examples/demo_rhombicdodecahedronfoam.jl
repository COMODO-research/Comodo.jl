using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

w = 1.0
n = (3,4,5)
E,V = rhombicdodecahedronfoam(w,n; merge=false, orientation=:allign)
F = element2faces(E)

## Visualize mesh
strokewidth = 2
strokecolor = :black
cmap = reverse(cgrad(:Spectral, length(E), categorical = true))

fig = Figure(size = (1200,1200))
# ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Rhombic dodecahedron")
ax1 = LScene(fig[1,1])
cc = Makie.Camera3D(ax1.scene, projectiontype = Makie.Orthographic)

hp1 = poly!(ax1, GeometryBasics.Mesh(V,F), color=:white,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)

fig
