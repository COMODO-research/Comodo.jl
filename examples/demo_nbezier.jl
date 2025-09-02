using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

# This demo shows the use of nbezier for Bezier curve interpolation. 

# Define control points
P = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
P[1 ] = GeometryBasics.Point{3, Float64}( 0.0, 0.0, 0.0)
P[2 ] = GeometryBasics.Point{3, Float64}( 1.0, 0.0, 0.0)
P[3 ] = GeometryBasics.Point{3, Float64}( 1.0, 1.0, 0.0)
P[4 ] = GeometryBasics.Point{3, Float64}( 1.0, 1.0, 1.0)

n = 25 # Number of points

V = nbezier(P,n) # Get Bezier fit points

# Visualization
GLMakie.closeall()

fig = Figure(size = (800,800))
ax = AxisGeom(fig[1, 1])

hp1 = scatter!(ax, P,markersize=25,color=:black)
scatter!(ax, V,markersize=15,color=:red)
hp2 = lines!(ax, V,linewidth=3,color=:red)

Legend(fig[1, 2],[hp1,hp2],["Control points","Bezier spline"])

fig
