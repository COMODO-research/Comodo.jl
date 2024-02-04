using Gibbon
using GLMakie
using GeometryBasics

r = 1.0 #radius
M1 = platonicsolid(1,r)
N1,VN1 = meshnormal(M1)

M2 = platonicsolid(2,r)
N2,VN2 = meshnormal(M2)

M3 = platonicsolid(3,r)
N3,VN3 = meshnormal(M3)

M4 = platonicsolid(4,r)
N4,VN4 = meshnormal(M4)

M5 = platonicsolid(5,r)
N5,VN5 = meshnormal(M5)

#Visualize mesh
fig = Figure(size = (1600,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Tetrahedron")
hp1 = poly!(ax1,M1, strokewidth=3,color=:red, shading=FastShading, overdraw=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Hexahedron (cube)")
hp2 = poly!(ax2,M2, strokewidth=3,color=:green, shading=FastShading, overdraw=false)

ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Octahedron")
hp3 = poly!(ax3,M3, strokewidth=3,color=:blue, shading=FastShading, overdraw=false)

ax4 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Icosahedron")
hp4 = poly!(ax4,M4, strokewidth=3,color=:yellow, shading=FastShading, overdraw=false)

ax5 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Dodecahedron")
hp5 = poly!(ax5,M5, strokewidth=3,color=:magenta, shading=FastShading, overdraw=false)

fig