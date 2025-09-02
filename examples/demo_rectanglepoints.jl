using Comodo
using Comodo.GLMakie

#=
This demo shows the use of `rectangle` function to create points on a rectangle
=#

# Input parameters 
w1 = 1.0 # width (x-direction)
h1 = 2.3 # height (y-direction)
dir1=:acw # Direction 

w2 = 1.0 # width (x-direction)
h2 = 2.3 # height (y-direction)
dir2=:cw # Direction 

w3 = 2.0 # width (x-direction)
h3 = 1.0 # height (y-direction)
dir3=:acw # Direction 
pointSpacing3 = 0.15 

w4 = 3.0 # width (x-direction)
h4 = 4.0 # height (y-direction)
dir4=:cw # Direction 
pointSpacing4 = 0.4 

# Create rectangle coordinates
V1 = rectanglepoints(w1, h1; dir=dir1) 
V2 = rectanglepoints(w2, h2; dir=dir2) 
V3 = rectanglepoints(w3, h3, pointSpacing3; dir=dir3) 
V4 = rectanglepoints(w4, h4, pointSpacing4; dir=dir4) 

# Visualization
GLMakie.closeall()

markersize1 = 10
markersize2 = 15

fig = Figure(size = (1600,1600))

ax1 = AxisGeom(fig[1, 1], title="Rectangle, w=$w1, h=$h1, dir=$dir1", azimuth=-pi/2, elevation=pi/2)
lines!(ax1, V1, linewidth=2.0, color=:red)
scatter!(ax1, V1, color=:black, markersize=markersize1)
scatter!(ax1, V1[1], color=:red, markersize=markersize2)

ax2 = AxisGeom(fig[1, 2], title="Rectangle, w=$w2, h=$h2, dir=$dir2", azimuth=-pi/2, elevation=pi/2)
lines!(ax2, V2, linewidth=2.0, color=:red)
scatter!(ax2, V2, color=:black, markersize=markersize1)
scatter!(ax2, V2[1], color=:red, markersize=markersize2)

ax3 = AxisGeom(fig[2, 1], title="Rectangle, w=$w3, h=$h3, pointSpacing=$pointSpacing3, dir=$dir3", azimuth=-pi/2, elevation=pi/2)
lines!(ax3, V3, linewidth=2.0, color=:red)
scatter!(ax3, V3, color=:black, markersize=markersize1)
scatter!(ax3, V3[1], color=:red, markersize=markersize2)

ax4 = AxisGeom(fig[2, 2], title="Rectangle, w=$w4, h=$h4, pointSpacing=$pointSpacing4, dir=$dir4", azimuth=-pi/2, elevation=pi/2)
lines!(ax4, V4, linewidth=2.0, color=:red)
scatter!(ax4, V4, color=:black, markersize=markersize1)
scatter!(ax4, V4[1], color=:red, markersize=markersize2)

fig