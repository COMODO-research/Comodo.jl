using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.DelaunayTriangulation

#=
This demo shows the use of the `erodeboundary!` function to obtain erode 
boundary triangles (i.e. remove them) based on a criteria. 
=#

GLMakie.closeall()

for testCase = 1:5
    if testCase == 1
        w = 1.0
        pointSpacing = 0.5
        V = squarepoints(w, pointSpacing)

        # Create unconstrained Delaunay triangulation to show difference with convex hull
        TR = triangulate(V; delete_ghosts=true)
        Fd = [TriangleFace{Int}(f) for f in each_solid_triangle(TR)] 
        Vd = get_points(TR)

        Fdc = Fd
        Vdc = Vd
        V = deepcopy(Vdc) # Copy so we keep original too for visualisation        

        α = pointSpacing
        method = :edgelength
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
        
        VT = (V,)
        R = ([1],)
        P = (pointSpacing)
        Fdc,Vdc,_ = regiontrimesh(VT,R,P)

        # Create unconstrained Delaunay triangulation to show difference with convex hull
        TR = triangulate(Vdc; delete_ghosts=true)
        Fd = [TriangleFace{Int}(f) for f in each_solid_triangle(TR)] 
        Vd = get_points(TR)

        Fdc = deepcopy(Fd)
        Vdc = deepcopy(Vd) # Copy so we keep original too for visualisation

        α = pointSpacing
        method = :edgelength
    elseif testCase == 3        
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

        VT = (V1, V2)
        R = ([1,2],)
        P = (pointSpacing)
        Fdc,Vdc,_ = regiontrimesh(VT,R,P)

        # Create unconstrained Delaunay triangulation to show difference with convex hull
        TR = triangulate(Vdc; delete_ghosts=true)
        Fd = [TriangleFace{Int}(f) for f in each_solid_triangle(TR)] 
        Vd = get_points(TR)

        Fdc = deepcopy(Fd)
        Vdc = deepcopy(Vd) # Copy so we keep original too for visualisation

        α = pointSpacing
        method = :edgelength
    elseif testCase == 4
        # Example point set (triangulated batman curve)
        n = 120
        V = batman(n; stepwise=true)
        pointSpacing = pointspacingmax(V)

        VT = (V,)
        R = ([1],)
        P = (pointSpacing)
        Fdc,Vdc,_ = regiontrimesh(VT,R,P)

        # Create unconstrained Delaunay triangulation to show difference with convex hull
        TR = triangulate(Vdc; delete_ghosts=true)
        Fd = [TriangleFace{Int}(f) for f in each_solid_triangle(TR)] 
        Vd = get_points(TR)

        Fdc = deepcopy(Fd)
        Vdc = deepcopy(Vd) # Copy so we keep original too for visualisation

        α = pointSpacing
        method = :circumcircle
    elseif testCase == 5 
        V = [Point{3,Float64}(0.0,0.0,0.0),
         Point{3,Float64}(1/3,0.0,0.0),
         Point{3,Float64}(2/3,0.0,0.0),
         Point{3,Float64}(1.0,0.0,0.0),
         Point{3,Float64}(1.0,1/3,0.0),
         Point{3,Float64}(2/3,1/3,0.0),
         Point{3,Float64}(2/3,2/3,0.0),
         Point{3,Float64}(1.0,2/3,0.0),
         Point{3,Float64}(1.0,1.0,0.0),
         Point{3,Float64}(2/3,1.0,0.0),
         Point{3,Float64}(1/3,1.0,0.0),
         Point{3,Float64}(0.0,1.0,0.0),
         Point{3,Float64}(0.0,2/3,0.0),
         Point{3,Float64}(0.0,1/3,0.0),
         ]

        V = subcurve(V,4; close_loop=true)

        α = 1.01*pointspacingmax(V; close_loop=true)

        Fd, Vd = delaunay2D(V) # Contruct Delaunay triangulation

        Fdc = deepcopy(Fd)
        Vdc = deepcopy(Vd) # Copy so we keep original too for visualisation

        method = :edgelength
    end

    # Compute alpha complex
    erodeboundary!(Fd,Vd,α; method=method)

    # Visualisation
    fig = Figure(size=(1200,1000))
    ax1 = AxisGeom(fig[1, 1], title="Delaunay triangulation", azimuth=-pi/2, elevation=pi/2)
    hp1 = meshplot!(ax1,Fdc,Vdc)    
    ax2 = AxisGeom(fig[1, 2], title="Eroded mesh, α=$α, method=$method", azimuth=-pi/2, elevation=pi/2)
    hp2 = meshplot!(ax2,Fd,Vd)

    # hp1 = edgeplot!(ax1,Eb,Vp; color=:blue, linewidth=3)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end