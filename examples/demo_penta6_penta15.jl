using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

testCase = 1

if testCase ==1 
    plateDim1 = [10.0,14.0]
    pointSpacing1 = 2.0

    orientation1 = :up
    F1,V1 = triplate(plateDim1,pointSpacing1; orientation=orientation1)
    t = 6
    n = 4
    direction=:positive
    E, V = extrudefaces(F1,V1; extent=t, direction=direction, num_steps=n)
    F = element2faces(E) # Triangular faces
end

E_penta15, V_penta15 = penta6_penta15(E,V)
F_penta15 = element2faces(E_penta15)

## Visualization
cmap = cgrad(:Spectral, 5, categorical = true)

strokewidth = 1 

fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Penta5 faces")
hp1 = poly!(ax1,GeometryBasics.Mesh(V,F[1]), color=:white, shading = FastShading, transparency=false, strokecolor=:black, strokewidth=1)
hp2 = poly!(ax1,GeometryBasics.Mesh(V,F[2]), color=:white, shading = FastShading, transparency=false, strokecolor=:black, strokewidth=1)
scatter!(ax1,V,color=:black,markersize=10)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Penta15 faces")
hp3 = poly!(ax2,GeometryBasics.Mesh(V_penta15,F_penta15[1]), color=:white, shading = FastShading, transparency=false, strokecolor=:black, strokewidth=1)
hp4 = poly!(ax2,GeometryBasics.Mesh(V_penta15,F_penta15[2]), color=:white, shading = FastShading, transparency=false, strokecolor=:black, strokewidth=1)
scatter!(ax2,V_penta15,color=:black,markersize=10)
# normalplot(ax2,F_penta15[2],V_penta15)

fig