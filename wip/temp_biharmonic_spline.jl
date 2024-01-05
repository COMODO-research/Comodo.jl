using Gibbon
using GLMakie

# Define raw data
x = range(0,9,9) # Interval definition
y = 5.0*cos.(x.^2 ./ 9.0) 
n = 50
xi = range(-0.5,9.5,n) # Interval definition

yi = interp_biharmonicSpline(x,y,xi; extrapolate_method="linear",pad_data="linear")

# Visualize 
fig1 = Figure()

ax = Axis(fig1[1, 1], aspect = DataAspect())

scatter!(ax, x,y,markersize=25,color=:black)
lines!(ax, x,y,linewidth=2,color=:black)

scatter!(ax, xi,yi,markersize=15,color=:red)
lines!(ax, xi,yi,linewidth=4,color=:red)

fig1
