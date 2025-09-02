using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

w = 1.0
n = (3,4,5)
E,V = rhombicdodecahedronfoam(w,n; merge=false, orientation=:align)
F = element2faces(E)

## Visualize mesh
GLMakie.closeall()

strokewidth = 2
strokecolor = :black
cmap = reverse(cgrad(:Spectral, length(E), categorical = true))

fig = Figure(size = (1200,1200))
ax1 = AxisGeom(fig[1, 1], title = "Rhombic dodecahedron")
hp1 = meshplot!(ax1, F, V)
display(fig)