using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This is a demonstration of the capabilities of the `tridisc` function which 
generates the faces `F` and vertices `V` for a triangulated disc (circle).
=#

# Define input parameters
r = 1.0 # Radius

# Create disc mesh 
n1 = 0
F1,V1 = tridisc(r,n1)

n2 = 2
F2,V2 = tridisc(r,n2)

n3 = 3
ngon3 = 5
method3 = :linear
F3,V3 = tridisc(r,n3; ngon=ngon3, method=method3)

n4 = 3
ngon4 = 6
method4 = :Loop
F4,V4 = tridisc(r,n4; ngon=ngon4, method=method4)

# Visualization
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Triangulated disc, n=$n1")
hp1 = poly!(ax1,GeometryBasics.Mesh(V1,F1), strokewidth = 1, color = :white, shading = FastShading, transparency = false)
# normalplot(ax1,F1,V1)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Triangulated disc, n=$n2")
hp2 = poly!(ax2,GeometryBasics.Mesh(V2,F2), strokewidth = 1, color = :white, shading = FastShading, transparency = false)

ax3 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Triangulated disc, n=$n3, ngon=$ngon3, method=$method3")
hp3 = poly!(ax3,GeometryBasics.Mesh(V3,F3), strokewidth = 1, color = :white, shading = FastShading, transparency = false)

ax4 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Triangulated disc, n=$n4, ngon=$ngon4, method=$method4")
hp4 = poly!(ax4,GeometryBasics.Mesh(V4,F4), strokewidth = 1, color = :white, shading = FastShading, transparency = false)

fig