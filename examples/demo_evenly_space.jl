using Comodo
using GeometryBasics
using GLMakie

# This demo shows how to evenly space a curve using spline interpolations.

testCase = 2

# Define input curve
if testCase == 1
    n = 5
    t = range(0,2π-2π/n,n)
    V = [GeometryBasics.Point{3, Float64}(3.0*cos(tt),3.0*sin(tt),0.0) for tt ∈ t]
elseif testCase == 2
    n = 4*7+4
    rFun(t) = sin(4*t)+2
    V = circlepoints(rFun,n)
end


# Evenly sample curve
n = 25; # Number of points for spline 
Vi1 = evenly_space(V, 0.3; niter = 10, close_loop = false) # Returns points and spline interpolation object

Vi2 = evenly_space(V, 0.3; niter = 10, close_loop = true) # Returns points and spline interpolation object


# Visualization

markersize1 = 10
markersize2 = 15
linewidth = 2

fig = Figure(size = (1200,500))

ax1 = Axis3(fig[1, 1],aspect = :data,azimuth=-pi/2,elevation=pi/2, title="Input")
hp1 = scatter!(ax1, V,markersize=markersize2,color=:red)
hp2 = lines!(ax1, V,linewidth=linewidth,color=:red)

ax2 = Axis3(fig[1, 2],aspect = :data,azimuth=-pi/2,elevation=pi/2, title="evenly spaced")
hp1 = scatter!(ax2, V,markersize=markersize2,color=:red)
hp2 = lines!(ax2, V,linewidth=linewidth,color=:red)
hp3 = scatter!(ax2, Vi1,markersize=markersize1,color=:black)
hp4 = lines!(ax2, Vi1,linewidth=linewidth,color=:black)

ax3 = Axis3(fig[1, 3],aspect = :data,azimuth=-pi/2,elevation=pi/2, title="evenly spaced, closed")
hp1 = scatter!(ax3, V,markersize=markersize2,color=:red)
hp2 = lines!(ax3, V,linewidth=linewidth,color=:red)
hp3 = scatter!(ax3, Vi2,markersize=markersize1,color=:black)
hp4 = lines!(ax3, Vi2,linewidth=linewidth,color=:black)

fig
