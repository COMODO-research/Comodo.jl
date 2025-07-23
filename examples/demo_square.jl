using Comodo
using Comodo.GLMakie

#=
This demo shows the use of `square` function to create mesh of a square
=#

# Input parameters 
w1 = 1.0 # width (x-direction)
orientation1=:up # Direction 

w2 = 6.0 # width (x-direction)
orientation2=:down # Direction 

# Create square coordinates
F1, V1 = square(w1; orientation=orientation1) 
F2, V2 = square(w2; orientation=orientation2) 

# Visualization
markersize1 = 10
markersize2 = 15
fig = Figure(size = (1600,1600))

ax1 = AxisGeom(fig[1, 1], title="square, w=$w1, orientation=$orientation1")
meshplot!(ax1, F1, V1; strokewidth=2.0)

ax2 = AxisGeom(fig[1, 2], title="square, w=$w2, orientation=$orientation2")
meshplot!(ax2, F2, V2; strokewidth=2.0)

fig