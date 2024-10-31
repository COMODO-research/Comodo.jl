using Comodo
using GLMakie
using GeometryBasics
using LinearAlgebra

r = 2.0 # Radius

n1 = 0 # Number of refinement steps from cube
Fn1,Vn1 = quadsphere(n1,1)

n2 = 1 # Number of refinement steps from cube
Fn2,Vn2 = quadsphere(n2,r)

n3 = 2 # Number of refinement steps from cube
Fn3,Vn3 = quadsphere(n3,r)

n4 = 3 # Number of refinement steps from cube
Fn4,Vn4 = quadsphere(n4,r)

## Visualization
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=$n1")
poly!(ax1,GeometryBasics.Mesh(Vn1,Fn1), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=$n2")
poly!(ax2,GeometryBasics.Mesh(Vn2,Fn2), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax3 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=$n3")
poly!(ax3,GeometryBasics.Mesh(Vn3,Fn3), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax4 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=$n4")
poly!(ax4,GeometryBasics.Mesh(Vn4,Fn4), strokewidth=3,color=:white,shading=FastShading,transparency=false)

fig