using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This is a demonstration of the capabilities of the `quaddisc` function which 
generates the faces `F` and vertices `V` for a quadrangulated disc (circle).
=#

# Define input parameters
r = 1.0 # Radius

n1 = 0
method1 = :linear
F1,V1 = quaddisc(r,n1; method=method1)

n2 = 1
method2 = :Catmull_Clark
F2,V2 = quaddisc(r,n2; method=method2)

n3 = 2
method3 = :Catmull_Clark
F3,V3 = quaddisc(r,n3; method=method3)

n4 = 3
method4 = :Catmull_Clark
F4,V4 = quaddisc(r,n4; method=method4, orientation=:down)

# Visualization
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Quandrangulated disc, n=$n1, method=$method1")
hp1 = poly!(ax1,GeometryBasics.Mesh(V1,F1), strokewidth=1,color=:white,shading=FastShading,transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Quandrangulated disc, n=$n2, method=$method2")
hp2 = poly!(ax2,GeometryBasics.Mesh(V2,F2), strokewidth=1,color=:white,shading=FastShading,transparency=false)

ax3 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Quandrangulated disc, n=$n3, method=$method3")
hp3 = poly!(ax3,GeometryBasics.Mesh(V3,F3), strokewidth=1,color=:white,shading=FastShading,transparency=false)

ax4 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Quandrangulated disc, n=$n4, method=$method4")
hp4 = poly!(ax4,GeometryBasics.Mesh(V4,F4), strokewidth=1,color=:white,shading=FastShading,transparency=false)

fig