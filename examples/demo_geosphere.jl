using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of the geosphere function. An unrefined sphere is an 
icosahedron. Through subdivision (see `subtri`) a refined geodesic dome is 
obtained. 
=#

r = 1.0 # radius
n1 = 0 # Number of refinement iterations
F1,V1 = geosphere(n1,r)

n2 = 1 # Number of refinement iterations
F2,V2 = geosphere(n2,r)

n3 = 2 # Number of refinement iterations
method3 = :linear
F3,V3 = geosphere(n3,r; method=method3)

n4 = 3 # Number of refinement iterations
method4 = :Loop
F4,V4 = geosphere(n4,r; method=method4)

#Visualize mesh
GLMakie.closeall()

fig = Figure(size = (1600,800))

ax1 = AxisGeom(fig[1, 1], title = "A geodesic sphere n=$n1")
hp1 = meshplot!(ax1, F1, V1)

ax2 = AxisGeom(fig[1, 2], title = "A geodesic sphere n=$n2")
hp2 = meshplot!(ax2, F2, V2)

ax3 = AxisGeom(fig[2, 1], title = "A geodesic sphere n=$n3, method=$method3")
hp3 = meshplot!(ax3, F3, V3)

ax4 = AxisGeom(fig[2, 2], title = "A geodesic sphere n=$n4, method=$method4")
hp4 = meshplot!(ax4, F4, V4)

fig