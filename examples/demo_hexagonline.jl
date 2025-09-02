using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 1.0
n = 25
V1 = hexagonline(r,n; type=:ufdf)
V2 = hexagonline(r,n; type=:zigzag)

## Visualize mesh
GLMakie.closeall()

markersize1 = 25
markersize2 = 15
linewidth1 = 2
linewidth2 = 1

fig = Figure(size = (1200,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "hexagonline, type=:ufdf")
hp1 = scatter!(ax1,V1,markersize=markersize2,color=:black)
hp2 = lines!(ax1,V1,linewidth=linewidth2,color=:black)

ax2 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "hexagonline, type=:zigzag")
hp1 = scatter!(ax2,V2,markersize=markersize2,color=:black)
hp2 = lines!(ax2,V2,linewidth=linewidth2,color=:black)

fig