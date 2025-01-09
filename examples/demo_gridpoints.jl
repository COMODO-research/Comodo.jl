using Comodo
using Comodo.GLMakie

#=
This demo shows the use of the gridpoints function. This function can be used to create
a grid of points. 
=#


# Define a set of points on a grid using a range, y and z are assumed equal to x when not provided
V1 = gridpoints(range(0.0,2.0*π,5)) 

# Define a set of points on a grid using a vector, y and z are assumed equal to x when not provided
V2 = gridpoints(collect(0:1:3)) 

# Define a set of points on a grid using a vector, z is assumed equal to x when not provided
V3 = gridpoints(collect(0:1:3),collect(0:1:5)) 

# Define a set of points on a grid using ranges
V4 = gridpoints(0.0:1.0:4.0,LinRange(-3.0,3.0,4),range(0,2.0*π,6)) 

# Visualize 
markersize =25

fig = Figure(size = (1200, 800))

ax1 = Axis3(fig[1, 1], aspect = :data)
scatter!(ax1, V1,markersize=markersize,color=[norm(v) for v in V1],colormap = :Spectral,strokewidth=2)

ax2 = Axis3(fig[1, 2], aspect = :data)
scatter!(ax2, V2,markersize=markersize,color=[norm(v) for v in V2],colormap = :Spectral,strokewidth=2)

ax3 = Axis3(fig[2, 1], aspect = :data)
scatter!(ax3, V3,markersize=markersize,color=[norm(v) for v in V3],colormap = :Spectral,strokewidth=2)

ax4 = Axis3(fig[2, 2], aspect = :data)
scatter!(ax4, V4,markersize=markersize,color=[norm(v) for v in V4],colormap = :Spectral,strokewidth=2)

fig