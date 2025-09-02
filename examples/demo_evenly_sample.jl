using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

# This demo shows how to evenly sample a curve using spline interpolations.

GLMakie.closeall()

for testCase = 1:3
    # Define input curve
    if testCase == 1
        n = 5
        t = range(0,2π-2π/n,n)
        V = [GeometryBasics.Point{3, Float64}(3.0*cos(tt),3.0*sin(tt),0.0) for tt ∈ t]

        # Evenly sample curve
        n = 50; # Number of points for spline 
        Vi1 = evenly_sample(V, n; spline_order=4, close_loop = false) # Returns points and spline interpolation object
        Vi2 = evenly_sample(V, n; close_loop = true) # Returns points and spline interpolation object
    elseif testCase == 2
        n = 4*7+4
        rFun(t) = sin(4*t)+2
        V = circlepoints(rFun,n)

        # Evenly sample curve
        n = 75; # Number of points for spline 
        Vi1 = evenly_sample(V, n; spline_order=4, close_loop = false, niter = 10) # Returns points and spline interpolation object
        Vi2 = evenly_sample(V, n; spline_order=4, close_loop = true, niter = 10) # Returns points and spline interpolation object
    elseif testCase == 3
        V = batman(25; stepwise = true, dir=:acw)

        # Evenly sample curve
        n = 50; # Number of points for spline 
        Vi1 = evenly_sample(V, n; niter = 10, spline_order=4, close_loop = false) # Returns points and spline interpolation object
        Vi2 = evenly_sample(V, n; niter = 10, spline_order=4, close_loop = true) # Returns points and spline interpolation object
    end

    # Visualization
    markersize1 = 10
    markersize2 = 15
    linewidth = 2

    fig = Figure(size = (1200,500))

    ax1 = AxisGeom(fig[1, 1], azimuth=-pi/2, elevation=pi/2, title="Input")
    hp1 = scatter!(ax1, V,markersize=markersize2,color=:red)
    hp2 = lines!(ax1, V,linewidth=linewidth,color=:red)

    ax2 = AxisGeom(fig[1, 2], azimuth=-pi/2,elevation=pi/2, title="evenly sampled")
    hp1 = scatter!(ax2, V,markersize=markersize2,color=:red)
    hp2 = lines!(ax2, V,linewidth=linewidth,color=:red)
    hp3 = scatter!(ax2, Vi1,markersize=markersize1,color=:black)
    hp4 = lines!(ax2, Vi1,linewidth=linewidth,color=:black)

    ax3 = AxisGeom(fig[1, 3], azimuth=-pi/2,elevation=pi/2, title="evenly sampled, closed")
    hp1 = scatter!(ax3, V,markersize=markersize2,color=:red)
    hp2 = lines!(ax3, V,linewidth=linewidth,color=:red)
    hp3 = scatter!(ax3, Vi2,markersize=markersize1,color=:black)
    hp4 = lines!(ax3, Vi2,linewidth=linewidth,color=:black)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end

