
using Gibbon
using GLMakie

testCase=1
if testCase==1
    M=icosahedron()
elseif testCase==2
    M=tetrahedron()
elseif testCase==3
    M=cube()
end

N,VN=meshnormal(M)

#Visualize mesh
GLMakie.activate!(inline=false) # To avoid plotting in plotpane as per: https://github.com/MakieOrg/Makie.jl/issues/2956
fig = Figure()

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A mesh with face normals visualized")
hp1=poly!(ax1,M, strokewidth=3,shading=true,color=:white, transparency=true, overdraw=false)
hpa=arrows!(ax1,VN,N)

fig
