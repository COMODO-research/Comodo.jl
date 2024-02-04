using Gibbon
using GLMakie
using GeometryBasics
using LinearAlgebra

r = 1.0 # Radius
n = 1 # Number of refinement steps from cube

Fn1,Vn1 = quadsphere(r,n)

Fn2,Vn2 = quadsphere(r,2)
Fn3,Vn3 = quadsphere(r,3)

## Visualization
fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=1")
poly!(ax1,GeometryBasics.Mesh(Vn1,Fn1), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=2")
poly!(ax2,GeometryBasics.Mesh(Vn2,Fn2), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=3")
hp2 = poly!(ax3,GeometryBasics.Mesh(Vn3,Fn3), strokewidth=3,color=:white,shading=FastShading,transparency=false)

fig