# using Gibbon
using GLMakie
using Interpolations

# Define raw data
x = 0:1:10 # Interval definition
y = 5.0*cos.(x.^2 ./ 9.0) 
z = sin.(x.^2 ./ 9.0)  

# Interpolation functions
# itp_linear = linear_interpolation(x, y, extrapolation_bc = Interpolations.Line())
# itp_cubic = cubic_spline_interpolation(x, y, extrapolation_bc = Interpolations.Line())

t = pushfirst!(cumsum(sqrt.(diff(x).^2 .+ diff(y).^2)),0.0)
ti = range(0.0,t[end],length(t)*5)

# Interpolate

x_new = lerp(t,x,ti)
y_new = lerp(t,y,ti)

# Visualisation
fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect())
scatter!(ax,x,y, color = :black, markersize = 25)
lines!(ax,x_new, y_new, linewidth = 3, color = :red,  linestyle = :solid)
scatter!(ax,x_new, y_new, color = :red, markersize = 15)

fig