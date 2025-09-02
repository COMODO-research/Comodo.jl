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
GLMakie.closeall()

fig = Figure(size=(800,800))

ax1 = AxisGeom(fig[1, 1], title = "Quandrangulated disc, n=$n1, method=$method1")
hp1 = meshplot!(ax1, F1, V1)

ax2 = AxisGeom(fig[1, 2], title = "Quandrangulated disc, n=$n2, method=$method2")
hp2 = meshplot!(ax2, F2, V2)

ax3 = AxisGeom(fig[2, 1], title = "Quandrangulated disc, n=$n3, method=$method3")
hp3 = meshplot!(ax3, F3, V3)

ax4 = AxisGeom(fig[2, 2], title = "Quandrangulated disc, n=$n4, method=$method4")
hp4 = meshplot!(ax4, F4, V4)

fig