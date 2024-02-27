using Comodo
using GLMakie
using BSplineKit

"""
This demo shows the use of the BSplineKit package for spline based curve interpolation. 
"""

# Define raw data
x = range(0,9,9) # Interval definition
y = 5.0*cos.(x.^2 ./ 9.0) 
n = 100
xi = range(-0.5,9.5,n) # Interval definition

# Interpolate/extrapolate using a Natural Cubic spline 
S = extrapolate(interpolate(x, y, BSplineOrder(4),Natural()),Smooth())

yi = S.(xi)

# Visualization
fig = Figure(size = (800,800))
ax = Axis(fig[1, 1], aspect = DataAspect())

hp1 = scatter!(ax, x,y,markersize=25,color=:black)
scatter!(ax, xi,yi,markersize=15,color=:red)
hp2 = lines!(ax, xi,yi,linewidth=3,color=:red)
Legend(fig[1, 2],[hp1,hp2],["Raw","Interpolated"])

fig
