using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.LinearAlgebra
 
GLMakie.closeall()

for testCase = 1:6
    if testCase == 1
        r = 2.0
        n = 8
        V1 = circlepoints(r,n; dir=:acw)  
        ind1 = collect(1:length(V1))
        N1 = fill(Vec{3,Float64}(0.0,0.0,1.0),length(ind1))

        V2 = V1
        ind2 = ind1
        N2 = fill(Vec{3,Float64}(0.0,0.0,-1.0),length(ind2))

        close_loop = true
        anglethreshold = 180.0

        showSurf = false # For visualisation purposes
    elseif testCase == 2 
        n = 12
        t = range(0.0,2*pi,n)
        V1 = [Point{3,Float64}(t,cos(t),0.0) for t in t]
        ind1 = collect(1:length(V1))
        N1 = fill(Vec{3,Float64}(0.0,0.0,1.0),length(ind1))
        
        V2 = V1
        ind2 = ind1
        N2 = fill(Vec{3,Float64}(0.0,0.0,-1.0),length(ind2))

        close_loop = false
        anglethreshold = 180.0

        showSurf = false # For visualisation purposes
    elseif testCase == 3
        n = 12
        V1 = collect(range(Point{3,Float64}(0.0,0.0,0.0),Point{3,Float64}(n,0.0,0.0),n))
        for i in 3:3:n
            V1[i] = Point{3,Float64}(V1[i][1],1.0,0.0)
        end
        ind1 = collect(1:length(V1))
        N1 = fill(Vec{3,Float64}(0.0,0.0,1.0),length(ind1))

        V2 = V1
        ind2 = ind1
        N2 = fill(Vec{3,Float64}(0.0,0.0,-1.0),length(ind2))

        close_loop = false
        anglethreshold = 140.0

        showSurf = false # For visualisation purposes
    elseif testCase == 4
        n = 13
        V1 = collect(range(Point{3,Float64}(0.0,0.0,0.0),Point{3,Float64}(n,0.0,0.0),n))
        for i in 2:2:n
            V1[i] = Point{3,Float64}(V1[i][1],1.0,0.0)
        end
        ind1 = collect(1:length(V1))
        N1 = fill(Vec{3,Float64}(0.0,0.0,1.0),length(ind1))

        V2 = V1
        ind2 = ind1
        N2 = fill(Vec{3,Float64}(0.0,0.0,-1.0),length(ind2))

        close_loop = false
        anglethreshold = 140.0

        showSurf = false # For visualisation purposes
    elseif testCase == 5
        r = 2.0
        n = 16
        V1 = circlepoints(r,n; dir=:acw)
        for i in 1:2:n  
            V1[i]/=2.0
        end
        ind1 = collect(1:length(V1))
        N1 = fill(Vec{3,Float64}(0.0,0.0,1.0),length(ind1))

        V2 = V1
        ind2 = ind1
        N2 = fill(Vec{3,Float64}(0.0,0.0,-1.0),length(ind2))

        close_loop = true
        anglethreshold = 90.0

        showSurf = false # For visualisation purposes
    elseif testCase == 6 
        n = 3
        r = 2.0
        F,V = geosphere(n,r)    
        VF = simplexcenter(F,V)
        L = [vf[3]<0.1 for vf in VF]        
        F1,V1,_ = remove_unused_vertices(F[L],V)    
        Eb = boundaryedges(F1)
        ind1 = edges2curve(Eb)
        ind1 = ind1[1:end-1]    
        NV = vertexnormal(F1,V1)
        N1 = NV[ind1]

        F2,V2,_ = remove_unused_vertices(F[.!L],V)
        Eb = boundaryedges(F2)
        ind2 = edges2curve(Eb)
        ind2 = ind2[1:end-1]    
        NV = vertexnormal(F2,V2)
        N2 = NV[ind2]

        close_loop = true
        anglethreshold = 120.0

        showSurf = true # For visualisation purposes
    end

    F1n = triangulateboundary(V1, ind1, N1, anglethreshold; deg = true, close_loop=close_loop)
    F2n = triangulateboundary(V2, ind2, N2, anglethreshold; deg = true, close_loop=close_loop)

    ## Visualization
    strokewidth = 1
    linewidth = 3

    fig = Figure(size=(1200,800))

    ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Triangulated boundary 1")
    hp1 = lines!(ax1,V1[ind1], linewidth=linewidth, color = :blue)
    hp2 = poly!(ax1,GeometryBasics.Mesh(V1,F1n), strokewidth = strokewidth, color = :red, shading = FastShading)
    if showSurf
        hp3 = poly!(ax1,GeometryBasics.Mesh(V1,F1), strokewidth = strokewidth, color = :white, shading = FastShading)
    end
    ax1 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Triangulated boundary 2")
    hp1 = lines!(ax1,V2[ind2], linewidth=linewidth, color = :blue)
    hp2 = poly!(ax1,GeometryBasics.Mesh(V2,F2n), strokewidth = strokewidth, color = :red, shading = FastShading)
    if showSurf
        hp3 = poly!(ax1,GeometryBasics.Mesh(V2,F2), strokewidth = strokewidth, color = :white, shading = FastShading)
    end

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end