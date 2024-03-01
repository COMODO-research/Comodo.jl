using Comodo
using GeometryBasics
using GLMakie

# This demo shows how to evenly sample a curve using spline interpolations.

# Define raw curve
t = range(0,2.0*pi,20)
V = [GeometryBasics.Point{3, Float64}(tt,3.0*sin(tt),0.0) for tt ∈ t]

# Evenly sample curve
n = 50; # Number of points for spline 
Vi, S = evenly_sample(V, n; niter = 10) # Returns points and spline interpolation object

# Visualization
fig = Figure(size = (1200,800))

ax1 = Axis3(fig[1, 1],aspect = :data)
hp1 = scatter!(ax1, V,markersize=25,color=:red)
hp2 = lines!(ax1, V,linewidth=3,color=:red)

ax2 = Axis3(fig[1, 2],aspect = :data)
hp1 = scatter!(ax2, Vi,markersize=20,color=:black)
hp2 = lines!(ax2, Vi,linewidth=3,color=:black)
hp2 = lines!(ax2, V,linewidth=2,color=:red)

fig
