using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.DelaunayTriangulation

#=
This demo shows the use of the `simplicialcomplex` function to obtain the 
simplicialcomplex for a set of points. This demo shows derivation of the 
simplicialcomplex as well as the unconstrained Delaunay Triangulation 
(which triangulates the convex hull) for comparison. 
=#

GLMakie.closeall()

for testCase = 1:5
    if testCase == 1
        w = 1.0
        pointSpacing = 0.5
        V = squarepoints(w, pointSpacing)
    elseif testCase == 2
        r1 = 10.0
        r2 = 20.0
        pointSpacing = 2.0
        t = range(pi, 0.0, spacing2numsteps(r1*pi,pointSpacing))
        V1 = [Point{3,Float64}(r1*cos(tt), r1*sin(tt), 0.0) for tt in t]
        V2 = [Point{3,Float64}(x, 0.0, 0.0) for x in range(r1, r2, spacing2numsteps(r2-r1,pointSpacing))]    
        t = range(0.0, pi, spacing2numsteps(r2*pi,pointSpacing))
        V3 = [Point{3,Float64}(r2*cos(tt), r2*sin(tt), 0.0) for tt in t]
        V4 = [Point{3,Float64}(x, 0.0, 0.0) for x in range(-r2, -r1, spacing2numsteps(r2-r1,pointSpacing))]
        V = V1
        append!(V,V2[2:end-1])
        append!(V,V3)
        append!(V,V4[2:end-1])
    elseif testCase == 3
        # Example point set (triangulated batman curve)
        n = 120
        V = batman(n; stepwise=true)
    elseif testCase == 4       
        r1m = 20.0
        r1Fun(t) = r1m + 3.0*sin(8.0*t)        
        r2m = 15.0
        r2Fun(t) = r2m + 3.0*sin(8.0*t)        
        pointSpacing = 1.0
        np1 = spacing2numsteps(2.0*pi*r1m,pointSpacing; close_loop=true)
        np2 = spacing2numsteps(2.0*pi*r2m,pointSpacing; close_loop=true)
        V1 = circlepoints(r1Fun, np1; dir=:acw)
        V1 = evenly_space(V1,pointSpacing)
        V2 = circlepoints(r2Fun, np2; dir=:acw)              
        V2 = evenly_space(V2,pointSpacing)

        V = [V1;V2]        
    end

    Fd, Vd = delaunay2D(V) # Contruct Delaunay triangulation
    
    # Visualisation
    fig = Figure(size=(1200,1000))
    ax1 = AxisGeom(fig[1, 1], title="Delaunay triangulation", azimuth=-pi/2, elevation=pi/2)
    hp1 = meshplot!(ax1,Fd,Vd; color=:lightgreen, strokewidth=0.5)
    scatter!(ax1, V, color=:black, markersize=10, depth_shift=-0.01f0)
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end