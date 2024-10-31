using Comodo
using GLMakie
using GeometryBasics
using Rotations



w = 1.0
n = (4,5,5)
E,V = kelvinfoam(w,n; merge=false)
F = element2faces(E)

## Visualize mesh
strokewidth = 2
strokecolor = :black
cmap = reverse(cgrad(:Spectral, length(E), categorical = true))

fig = Figure(size = (1200,1200))
# ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Rhombic dodecahedron")
ax1 = LScene(fig[1,1]); cc = Makie.Camera3D(ax1.scene, projectiontype = Makie.Perspective)

hp1 = poly!(ax1, GeometryBasics.Mesh(V,F[1]), color=:white,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)
hp2 = poly!(ax1, GeometryBasics.Mesh(V,F[2]), color=:white,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)

display(fig)

cc.near[] = 1f-3
cc.far[] = 100
