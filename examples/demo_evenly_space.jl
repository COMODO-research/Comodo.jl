using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

# This demo shows how to evenly space a curve using spline interpolations.

testCase = 1

# Define input curve
if testCase == 1
    n = 5
    t = range(0,2π-2π/n,n)
    V = [GeometryBasics.Point{3, Float64}(3.0*cos(tt),3.0*sin(tt),0.0) for tt ∈ t]

    must_points = [1,2,3,4,5]

    # Evenly sample curve
    Vi1 = evenly_space(V, 0.65; close_loop = false, spline_order = 4) # Returns points and spline interpolation object
    Vi2 = evenly_space(V, 0.65; close_loop = true, spline_order = 4,must_points = must_points) # Returns points and spline interpolation object
elseif testCase == 2
    n = 4*7+4
    rFun(t) = sin(4*t)+2
    V = circlepoints(rFun,n)

    must_points = [3,7,11,15,19,23,27,31]
    
    # Evenly sample curve
    Vi1 = evenly_space(V, 0.5; close_loop = false, niter=10, spline_order = 4) # Returns points and spline interpolation object
    Vi2 = evenly_space(V, 0.3; close_loop = true, niter=10, spline_order = 4,must_points = must_points) # Returns points and spline interpolation object
elseif testCase == 3
    V = batman(25; symmetric = true, dir=:acw)

    must_points = [7,10,13,15,18,21]

    # Evenly sample curve    
    Vi1 = evenly_space(V, 0.15; close_loop = false, niter=10, spline_order = 4) # Returns points and spline interpolation object
    Vi2 = evenly_space(V, 0.1; close_loop = true, niter=10, spline_order = 4, must_points = must_points) # Returns points and spline interpolation object
end


# Visualization
markersize1 = 10
markersize2 = 15
markersize3 = 20

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
if !isnothing(must_points)
    hp5 = scatter!(ax3, V[must_points],markersize=markersize3,color=:green)
end
hp2 = lines!(ax3, V,linewidth=linewidth,color=:red)
hp3 = scatter!(ax3, Vi2,markersize=markersize1,color=:black)
hp4 = lines!(ax3, Vi2,linewidth=linewidth,color=:black)

fig
