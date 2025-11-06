using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 2.0
Fb,Vb = geosphere(1,r)
E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb)   


tetFaces = element2faces(E)
tetFacesUnique, indReverseFaces = gunique(tetFaces; return_unique=Val(true), return_inverse=Val(true), sort_entries=true)

E_tet15, V_tet15 = tet4_tet15(E,V)
F_tet15 = element2faces(E_tet15)

# Visualization
GLMakie.closeall()

cmap = cgrad(:Spectral, 5, categorical = true)

F = element2faces(E) # 6-noded triangle faces

fig = Figure(size=(800,800))

ax1 = AxisGeom(fig[1, 1], title = "Tet4 faces (3-noded)")
hp1 = meshplot!(ax1, F, V)
scatter!(ax1, V, color=:black, markersize=15)
# normalplot(ax1, F, V)

ax2 = AxisGeom(fig[1, 2], title = "Tet15 faces (6-noded)")
hp2 = meshplot!(ax2, F_tet15, V_tet15)
scatter!(ax2, V_tet15, color=:black, markersize=15)
# normalplot(ax2, F_tet10, V_tet10)
fig