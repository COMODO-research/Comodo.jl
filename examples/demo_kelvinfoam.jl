using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

w = 1.0
n = (3,4,5)
E,V = kelvinfoam(w,n; merge=false)
F = element2faces(E)

## Visualize mesh
GLMakie.closeall()

elementLabels = collect(1:length(E))
C1 = repeat(elementLabels,inner=8)
C2 = repeat(elementLabels,inner=6)

Fs1, Vs1 = separate_vertices(F[1], V)
Cs1_Vs1 = simplex2vertexdata(Fs1, C1)
Fs2, Vs2 = separate_vertices(F[2], V)
Cs2_Vs2 = simplex2vertexdata(Fs2, C2)

strokewidth = 2
strokecolor = :black
cmap = reverse(cgrad(:Spectral, length(E), categorical = true))

fig = Figure(size = (1200,1200))
ax1 = AxisGeom(fig[1, 1], title = "Rhombic dodecahedron")
hp1 = meshplot!(ax1, Fs1, Vs1, color=Cs1_Vs1, colormap=cmap)
hp2 = meshplot!(ax1, Fs2, Vs2, color=Cs2_Vs2, colormap=cmap)
display(fig)