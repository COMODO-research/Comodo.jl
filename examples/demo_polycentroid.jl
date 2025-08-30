using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.Statistics

GLMakie.closeall()

for testCase = 1:2
    if testCase == 1
        V = circlepoints(1.0,4)
        V2 = collect(range(V[end],V[1],10))
        append!(V,V2[2:end-1])
    elseif testCase == 2
        ## Defining a set of points on a circle using scalar radius
        r = 1.0 # Radius
        n = 40 # Number of points
        t = range(0.0,pi,11)
        V = r.*[Point{3,Float64}(cos(t),sin(t),0.0) for t in t]
        t = range(pi,2.0*pi,42)
        V2 = r.*[Point{3,Float64}(cos(t),sin(t),0.0) for t in t]
        append!(V,V2[2:end-1])
        Q = RotXYZ(0.0,-0.25*π,0.25*π) # Define a rotation tensor using Euler angles
        V = [Point{3, Float64}(Q*v) for v ∈ V] # Rotate the coordinates
    end

    C_true = Point{3, Float64}(2.1,3.4,-5.7) #randn(eltype(V))
    V .+= C_true # Shift coordinates to desired centre 

    # Compute polygon centroid
    C = polycentroid(V)

    ## Visualization
    fig = Figure(size=(900,900))
    ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")
    h0 = lines!(ax1, V, linewidth = 5, color = :black)
    h1 = scatter!(ax1, V, markersize=25, color = :black)
    h2 = scatter!(ax1, C_true, markersize=35, color = :green)
    h3 = scatter!(ax1, mean(V), markersize=25, color = :red)
    h4 = scatter!(ax1, C, markersize=25, color = :blue)
    Legend(fig[1, 2], [h1, h2, h3, h4],["Polygon points", "True centroid", "Average of points", "Derived centroid"])
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end