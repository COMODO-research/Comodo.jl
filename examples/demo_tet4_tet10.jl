using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 1.0
Fb,Vb = geosphere(2,r)
E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb)   




E_tet10, V_tet10 = tet4_tet10(E,V)
F_tet10 = element2faces(E_tet10)

## Visualization
cmap = cgrad(:Spectral, 5, categorical = true)

F = element2faces(E) # Triangular faces

strokewidth = 1 

fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Tet4 faces (tri3)")
hp1 = poly!(ax1,GeometryBasics.Mesh(V,F), color=:white, shading = FastShading, transparency=false, strokecolor=:black, strokewidth=1)
scatter!(ax1,V,color=:black,markersize=15)
# normalplot(ax1,F,V)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Tet10 faces (tri6)")
# hp1 = poly!(ax2,GeometryBasics.Mesh(V,Fb), color=:white, shading = FastShading, transparency=true, strokecolor=:black, strokewidth=1)
hp2 = poly!(ax2,GeometryBasics.Mesh(V_tet10,F_tet10), color=:white, shading = FastShading, transparency=false, strokecolor=:black, strokewidth=1)
# scatter!(ax1,V,color=:black,markersize=10)
# normalplot(ax2,F_tet10,V_tet10)
# scatter!(ax2,V,color=:black,markersize=15)
# scatter!(ax2,Vn,color=:red,markersize=15)
fig