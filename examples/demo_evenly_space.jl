using Comodo
using GeometryBasics
using GLMakie

# This demo shows how to evenly sample a curve using spline interpolations.

# Define raw curve
n = 5
t = range(0,2π-2π/n,n)
# t = range(0,2π,n)

V = [GeometryBasics.Point{3, Float64}(3.0*cos(tt),3.0*sin(tt),0.0) for tt ∈ t]

# Evenly sample curve
n = 50; # Number of points for spline 
Vi1 = evenly_space(V, 0.3; niter = 10, close_loop = false) # Returns points and spline interpolation object

Vi2 = evenly_space(V, 0.3; niter = 10, close_loop = true) # Returns points and spline interpolation object


# Visualization
fig = Figure(size = (1200,500))

ax1 = Axis3(fig[1, 1],aspect = :data,azimuth=-pi/2,elevation=pi/2, title="Input")
hp1 = scatter!(ax1, V,markersize=25,color=:red)
hp2 = lines!(ax1, V,linewidth=3,color=:red)

ax2 = Axis3(fig[1, 2],aspect = :data,azimuth=-pi/2,elevation=pi/2, title="evenly spaced")
hp1 = scatter!(ax2, V,markersize=25,color=:red)
hp2 = lines!(ax2, V,linewidth=3,color=:red)
hp3 = scatter!(ax2, Vi1,markersize=15,color=:black)
hp4 = lines!(ax2, Vi1,linewidth=3,color=:black)

ax3 = Axis3(fig[1, 3],aspect = :data,azimuth=-pi/2,elevation=pi/2, title="evenly spaced, closed")
hp1 = scatter!(ax3, V,markersize=25,color=:red)
hp2 = lines!(ax3, V,linewidth=3,color=:red)
hp3 = scatter!(ax3, Vi2,markersize=15,color=:black)
hp4 = lines!(ax3, Vi2,linewidth=3,color=:black)

fig
