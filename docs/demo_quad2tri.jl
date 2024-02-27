using Comodo
using GLMakie
using GeometryBasics
using LinearAlgebra

# Example data 
r = 1.0 # Radius
n = 3 # Number of refinement steps from cube
F,V = quadsphere(n,r)

convert_method = "angle"
Ft = quad2tri(F,V; convert_method = convert_method)

## Visualization
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original")
poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Forward slash converted")
poly!(ax2,GeometryBasics.Mesh(V,quad2tri(F,V; convert_method = "forward")), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax3 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Back slash converted")
poly!(ax3,GeometryBasics.Mesh(V,quad2tri(F,V; convert_method = "backward")), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax4 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Angle based conversion")
poly!(ax4,GeometryBasics.Mesh(V,quad2tri(F,V; convert_method = "angle")), strokewidth=3,color=:white,shading=FastShading,transparency=false)

fig