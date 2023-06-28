
using Gibbon
using GLMakie
using GeometryBasics

r=1.0 #radius
M1=platonicsolid(1,r)
N1,VN1=meshnormal(M1)

M2=platonicsolid(2,r)
N2,VN2=meshnormal(M2)

M3=platonicsolid(3,r)
N3,VN3=meshnormal(M3)

M4=platonicsolid(4,r)
N4,VN4=meshnormal(M4)

M5=platonicsolid(5,r)
N5,VN5=meshnormal(M5)

#Visualize mesh
GLMakie.activate!(inline=false) # To avoid plotting in plotpane as per: https://github.com/MakieOrg/Makie.jl/issues/2956
fig = Figure()

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Tetrahedron")
hp1=poly!(ax1,M1, strokewidth=3,shading=true,color=:red, transparency=true, overdraw=false)

ax2=Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Hexahedron (cube)")
hp2=poly!(ax2,M2, strokewidth=3,shading=true,color=:green, transparency=true, overdraw=false)

ax3=Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Octahedron")
hp3=poly!(ax3,M3, strokewidth=3,shading=true,color=:blue, transparency=true, overdraw=false)

ax4=Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Icosahedron")
hp4=poly!(ax4,M4, strokewidth=3,shading=true,color=:yellow, transparency=true, overdraw=false)

ax5=Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Dodecahedron")
hp5=poly!(ax5,M5, strokewidth=3,shading=true,color=:magenta, transparency=true, overdraw=false)

ax6=Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "All")
hp61=poly!(ax6,M1, strokewidth=3,shading=true,color=:red, transparency=true, overdraw=false)
hp62=poly!(ax6,M2, strokewidth=3,shading=true,color=:green, transparency=true, overdraw=false)
hp63=poly!(ax6,M3, strokewidth=3,shading=true,color=:blue, transparency=true, overdraw=false)
hp64=poly!(ax6,M4, strokewidth=3,shading=true,color=:yellow, transparency=true, overdraw=false)
hp65=poly!(ax6,M5, strokewidth=3,shading=true,color=:magenta, transparency=true, overdraw=false)

fig