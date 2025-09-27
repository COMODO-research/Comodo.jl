using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics

w = 1.0
n = (2,2,3)
E,V = rhombicdodecahedronfoam(w,n; merge=false, orientation=:align)
F = element2faces(E)


Eh,Vh = rhombicdodecahedron2hex(E,V)
Eh,Vh = separate_vertices(Eh, Vh; scaleFactor=1.0)
Fh = element2faces(Eh)

## Visualize mesh
GLMakie.closeall()

strokewidth = 2
strokecolor = :black
cmap = reverse(cgrad(:Spectral, length(E), categorical = true))

fig = Figure(size = (1200,1200))
ax1 = AxisGeom(fig[1, 1], title = "Rhombic dodecahedron")
hp1 = meshplot!(ax1, F, V,  color=(:white, 0.5), transparency=true)
# hpa = normalplot(ax1, F, V)

ax2 = AxisGeom(fig[1, 2], title = "Rhombic dodecahedron")
hp2 = meshplot!(ax2, Fh, Vh,  color=(:white, 0.5), transparency=true)
# hpa = normalplot(ax2, Fh, Vh)

display(fig);