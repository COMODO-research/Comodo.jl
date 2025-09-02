using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 1.0
Fb,Vb = geosphere(2,r)
E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb)   

E_tet10, V_tet10 = tet4_tet10(E,V)
F_tet10 = element2faces(E_tet10)


# Visualization
GLMakie.closeall()

cmap = cgrad(:Spectral, 5, categorical = true)

F = element2faces(E) # 6-noded triangle faces

fig = Figure(size=(800,800))

ax1 = AxisGeom(fig[1, 1], title = "Tet4 faces (tri3)")
hp1 = meshplot!(ax1, F, V)
scatter!(ax1, V, color=:black, markersize=15)
# normalplot(ax1, F, V)

ax2 = AxisGeom(fig[1, 2], title = "Tet10 faces (tri6)")
hp2 = meshplot!(ax2, F_tet10, V_tet10)
scatter!(ax2, V_tet10, color=:black, markersize=10)
# normalplot(ax2, F_tet10, V_tet10)
fig