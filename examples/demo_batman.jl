using Comodo
using Comodo.GLMakie

#=
This demo shows the use of the `batman` function to create a Batman logo curve. 
The curve is useful for testing surface meshing algorithms since it contains
sharp transitions and pointy features. 
=#

n = 100
V1 = batman(n; stepwise = true, dir=:acw)
V2 = batman(n; stepwise = false, dir=:cw)

# Visualisation
GLMakie.closeall()

fig = Figure(size=(1000,500))

ax1 = AxisGeom(fig[1, 1]; title="stepwise=true, approximately n points, anti-clockwise", azimuth=-pi/2, elevation=pi/2)
hp1 = lines!(ax1, V1,linewidth=3,color=:blue)
hp2 = scatter!(ax1, V1,markersize=8,color=:red)
hp2 = scatter!(ax1, V1[1],markersize=15,color=:yellow)
hp2 = scatter!(ax1, V1[2],markersize=15,color=:orange)

ax2 = AxisGeom(fig[1, 2]; title = "stepwise=false, n points exactly, clockwise", azimuth=-pi/2, elevation=pi/2)
hp1 = lines!(ax2, V2,linewidth=3,color=:blue)
hp2 = scatter!(ax2, V2,markersize=8,color=:red)
hp2 = scatter!(ax2, V2[1],markersize=15,color=:yellow)
hp2 = scatter!(ax2, V2[2],markersize=15,color=:orange)

fig