
using Gibbon
using GLMakie
using GeometryBasics

function testFunc(a=1)
    return a
end

M1=tetrahedron()
M2=octahedron()
M3=cube()
M4=icosahedron()
M5=dodecahedron()

#Visualize mesh
GLMakie.activate!(inline=false) # To avoid plotting in plotpane as per: https://github.com/MakieOrg/Makie.jl/issues/2956
fig = Figure()

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Tetrahedron")
hp1=poly!(ax1,M1, strokewidth=3,shading=true,color=:red, transparency=true, overdraw=false)

ax2=Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Octahedron")
hp2=poly!(ax2,M2, strokewidth=3,shading=true,color=:green, transparency=true, overdraw=false)

ax3=Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Hexahedron (cube)")
hp3=poly!(ax3,M3, strokewidth=3,shading=true,color=:blue, transparency=true, overdraw=false)

ax4=Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Icosahedron")
hp4=poly!(ax4,M4, strokewidth=3,shading=true,color=:yellow, transparency=true, overdraw=false)

ax5=Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Dodecahedron")
hp5=poly!(ax5,M5, strokewidth=3,shading=true,color=:magenta, transparency=true, overdraw=false)

fig