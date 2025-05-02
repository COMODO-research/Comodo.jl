using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Statistics

r = 10.0/pi # Radius

n1 = 0 # Number of refinement steps from cube
Fn1,Vn1 = subquadsphere(n1,1)

n2 = 1 # Number of refinement steps from cube
Fn2,Vn2 = subquadsphere(n2,r)

n3 = 2 # Number of refinement steps from cube
Fn3,Vn3 = subquadsphere(n3,r)

n4 = 3 # Number of refinement steps from cube
Fn4,Vn4 = subquadsphere(n4,r)

pointSpacing = 0.25

## Visualization
strokewidth = 2

fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=$n1")
poly!(ax1,GeometryBasics.Mesh(Vn1,Fn1), strokewidth=strokewidth,color=:white,shading=FastShading,transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=$n2")
poly!(ax2,GeometryBasics.Mesh(Vn2,Fn2), strokewidth=strokewidth,color=:white,shading=FastShading,transparency=false)

ax3 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=$n3")
poly!(ax3,GeometryBasics.Mesh(Vn3,Fn3), strokewidth=strokewidth,color=:white,shading=FastShading,transparency=false)

ax4 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Refined n=$n4")
poly!(ax4,GeometryBasics.Mesh(Vn4,Fn4), strokewidth=strokewidth,color=:white,shading=FastShading,transparency=false)

fig