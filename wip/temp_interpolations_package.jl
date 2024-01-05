# using Gibbon
using GLMakie
using Interpolations

# Based on: 
# https://github.com/JuliaMath/Interpolations.jl/blob/master/docs/src/convenience-construction.md#example-with-plotsjl

# Define raw data
x = (0:1.0:10) # Interval definition
y = cos.(x.^2 ./ 9.0) 

# Interpolation functions
itp_linear = linear_interpolation(x, y, extrapolation_bc = Interpolations.Line())
itp_cubic = cubic_spline_interpolation(x, y, extrapolation_bc = Interpolations.Line())

# Interpolate
x_new = collect(-1:0.1:11) # Denser interval
y_new_lin = itp_linear(x_new)
y_new_cub = itp_cubic(x_new)

# Visualisation
fig = Figure()
ax =Axis(fig[1, 1])
scatter!(ax,x,y, color = :black, markersize = 25)
lines!(ax,x_new, y_new_lin, linewidth = 3, color = :red,  linestyle = :solid)
lines!(ax,x_new, y_new_cub, linewidth = 3, color = :blue,  linestyle = :solid)
fig