using Comodo
using GLMakie
using GeometryBasics

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

n3 = 3 # Number of refinement iterations
method3 = :linear
F3,V3 = geosphere(n3,r; method=method3)

n4 = 3 # Number of refinement iterations
method4 = :Loop
F4,V4 = geosphere(n4,r; method=method4)


#Visualize mesh
lineWidth = 1

fig = Figure(size=(1600,800))

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A geodesic sphere n=$n1")
hp1=poly!(ax1,V1,F1, strokewidth=lineWidth,color=:white, shading = FastShading)

ax2=Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A geodesic sphere n=$n2")
hp2=poly!(ax2,V2,F2, strokewidth=lineWidth,color=:white, shading = FastShading)

ax3=Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A geodesic sphere n=$n3, method=$method3")
hp3=poly!(ax3,V3,F3, strokewidth=lineWidth,color=:white, shading = FastShading)

ax4=Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A geodesic sphere n=$n4, method=$method4")
hp4=poly!(ax4,V4,F4, strokewidth=lineWidth,color=:white, shading = FastShading)

fig