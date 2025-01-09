using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Random

#=
In this demo biharmonic interpolation is used for 3D scattered data 
interpolation. First a set of random 3D points and "value data" is defined. Next
a grid of points are defined onto which this value data is interpolated. 
=#

Random.seed!(1) # Set seed so demo performs the same each time

# Define raw data
m = 150 # Number of input data points
V = Vector{GeometryBasics.Point{3, Float64}}(undef,m) # Initialize point vector
C = Vector{Float64}(undef,m) # Initialize value/color vector
for q = 1:m # Loop over all points
    V[q] = 2.0*pi*rand(3) # Assign a randon point from 0-2*pi 
    C[q] = sin(2*sqrt(sum(V[q].^2))) # An interesting function with distance from origin
end

# Define points for interpolation 
n = 25
Vi = gridpoints(range(0.0,2.0*pi,n),range(0.0,2.0*pi,n),range(0,2.0*pi,5)) # Define a grid of points 

# Interpolate the data onto the grid
Ci = interp_biharmonic(V,C,Vi)

# Visualize 
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1], aspect = :data)

hs1 = scatter!(ax, V,markersize=35,color=C,colormap = :Spectral,strokewidth=2)
hs2 = scatter!(ax, Vi,markersize=15,color=Ci,colormap = :Spectral)
Colorbar(fig[1, 2],hs1,label = "Value data")
Legend(fig[1, 3],[hs1,hs2],["Raw","Interpolated"])
fig
