using Comodo
using Comodo.GLMakie

#=
This demo shows the use of `rectangle` function to create mesh of a rectangle
=#

# Input parameters 
w1 = 1.0 # width (x-direction)
h1 = 2.3 # height (y-direction)
orientation1=:up # Direction 

w2 = 6.0 # width (x-direction)
h2 = 4.0 # height (y-direction)
orientation2=:down # Direction 

# Create rectangle coordinates
F1, V1 = rectangle(w1, h1; orientation=orientation1) 
F2, V2 = rectangle(w2, h2; orientation=orientation2) 

# Visualization
GLMakie.closeall()

markersize1 = 10
markersize2 = 15

fig = Figure(size = (1600,1600))
ax1 = AxisGeom(fig[1, 1], title="rectangle, w=$w1, h=$h1, orientation=$orientation1")
meshplot!(ax1, F1, V1; strokewidth=2.0)
ax2 = AxisGeom(fig[1, 2], title="rectangle, w=$w2, h=$h2, orientation=$orientation2")
meshplot!(ax2, F2, V2; strokewidth=2.0)
fig