
using Gibbon
using GLMakie
using GeometryBasics

"""
This demo shows the use of the geoSphere function. An unrefined sphere is an 
icosahedron. Through subdivision (see `subtri`) a refined geodesic dome is 
obtained. 
"""

r = 1.0 # radius
n1 = 0 # Number of refinement iterations
F1,V1 = geoSphere(n1,r)

n2 = 1 # Number of refinement iterations
F2,V2 = geoSphere(n2,r)

n3 = 3 # Number of refinement iterations
F3,V3 = geoSphere(n3,r)

#Visualize mesh
fig = Figure(size=(1600,800))

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A geodesic sphere n=0")
hp1=poly!(ax1,V1,F1, strokewidth=2,color=:white, shading = FastShading)

ax2=Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A geodesic sphere n=1")
hp2=poly!(ax2,V2,F2, strokewidth=2,color=:white, shading = FastShading)

ax3=Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A geodesic sphere n=3")
hp3=poly!(ax3,V3,F3, strokewidth=2,color=:white, shading = FastShading)

fig