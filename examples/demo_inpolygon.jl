using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.DelaunayTriangulation

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1    
        w = 0.1
        V = [Point{3,Float64}(-w, 0.0, 0.0), Point{3,Float64}(-1-w, 1.0, 0.0), 
            Point{3,Float64}(1+w, 1.0, 0.0), Point{3,Float64}(w, 0.0, 0.0),
            Point{3,Float64}(1+w, -1.0, 0.0), Point{3,Float64}(-1-w, -1.0, 0.0)]#circlepoints(1.0,4; dir=:acw)
        
        np=7
        xRange = range(-1.0-w,1.0+w,np)
        yRange = range(-1.0,1.0,np)
        Vq = gridpoints(xRange, yRange,[0.0])   
    elseif testCase == 2
        w = 0.1
        V = [Point{3,Float64}(-w, 0.0, 0.0), Point{3,Float64}(-1-w, 1.0, 0.0), 
            Point{3,Float64}(1+w, 1.0, 0.0), Point{3,Float64}(w, 0.0, 0.0),
            Point{3,Float64}(1+w, -1.0, 0.0), Point{3,Float64}(-1-w, -1.0, 0.0)]#circlepoints(1.0,4; dir=:acw)
        n = 3
        V = subcurve(V,n; close_loop = true) 
        np = 7  
        Vq = 0.5.* randn(Point{2,Float64},500)    
    elseif testCase == 3
        n = 100
        r = 1.5
        rFun(t) = r + r/2.0 *sin(3*t)
        V = circlepoints(rFun,n; dir=:acw)
        pointSpacing = r/5.0
        V = evenly_space(V,pointSpacing; close_loop=true)

        # Random points
        Vq = randn(Point{2,Float64},500)            
    end

    # Some co-linear on edge points 
    Vq_add = [(0.25*V[i] + 0.75*V[i+1]) for i in 1:length(V)-1]
    append!(Vq,Vq_add)
    Vq_add = [(0.75*V[i] + 0.25*V[i+1]) for i in 1:length(V)-1]
    append!(Vq,Vq_add)

    # Some co-linear off edge points
    Vq_add = [(V[i] + 1.25*(V[i+1]-V[i])) for i in 1:length(V)-1]
    append!(Vq,Vq_add)
    Vq_add = [(V[i] - 0.25*(V[i+1]-V[i])) for i in 1:length(V)-1]
    append!(Vq,Vq_add)

    # Add offset points (same y-direction)
    Vq_add = [v.+Point{3,Float64}(-0.5, 0.0, 0.0) for v in V]
    append!(Vq,Vq_add)

    #Add copy of polygon points (fully same coordinates)
    append!(Vq,deepcopy(V))

    F = [inpolygon(p,V) for p in Vq]

    ## Visualization
    fig = Figure(size=(1600,800))
    # ax1 = LScene(fig[1,1]); cc = Makie.Camera3D(ax1.scene, projectiontype = Makie.Orthographic)
    ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Inpolygon testing")

    lines!(ax1,V,linewidth = 2, color = :blue)
    scatter!(ax1,V,markersize=18,color = :black)
    scatter!(ax1,Vq[F.==1],markersize=16,color = :green)
    scatter!(ax1,Vq[F.==0],markersize=12,color = :orange)
    scatter!(ax1,Vq[F.==-1],markersize=10,color = :red)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end