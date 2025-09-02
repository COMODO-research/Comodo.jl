using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO

GLMakie.closeall()

function plot_joingeom(F,V,C, testCase)
    c = Makie.Categorical(Makie.Reverse(:Spectral))

    Fn,Vn = separate_vertices(F,V)
    Cn = simplex2vertexdata(Fn,C,Vn)

    # Visualization
    fig = Figure(size=(1200,1200))

    ax1 = AxisGeom(fig[1, 1], title = "A multi-object mesh")
    hp1 = meshplot!(ax1, Fn, Vn, color=Cn, colormap=c)

    Legend(fig[1, 2], [hp1], ["Joined mesh"])
    Colorbar(fig[1, 3], hp1, label = "Labelling")

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end

for testCase = 1:3
    if testCase == 1 # Single face/point set 
        V = Vector{Point{3, Float64}}(undef,3)
        V[1 ] = Point{3, Float64}(0.0, 0.0, 0.0)
        V[2 ] = Point{3, Float64}(1.0, 0.0, 0.0)
        V[3 ] = Point{3, Float64}(1.0, 1.0, 0.0)
        F = Vector{TriangleFace{Int}}(undef,1)
        F[1 ] = TriangleFace{Int}(1,2,3)
        F,V = subtri(F,V,2)            

        F, V, C = joingeom(F, V)

        plot_joingeom(F,V,C, testCase)
    elseif testCase==2
        r = (3.0, 2.0, 1.0)
        F1,V1 = geosphere(2,r[1])            
        F2,V2 = geosphere(2,r[2])    
        V2 = [v.+Point{3,Float64}(r[1]+r[2],0.0,0.0) for v in V2]
        F3,V3 = geosphere(2,r[3])    
        V3 = [v.+Point{3,Float64}(-r[1]-r[3],0.0,0.0) for v in V3]

        F, V, C = joingeom(F1, V1, F2, V2, F3, V3)

        plot_joingeom(F,V,C, testCase)
    elseif testCase==3        
        F1,V1 = geosphere(2,1.0)            
        
        boxDim = [2.5,3.1,4] # Dimensions for the box in each direction
        pointSpacing = 0.5
        F2,V2,C2 = tribox(boxDim,pointSpacing)
        V2 = [v.+Point{3,Float64}(1.0+boxDim[1]/2+1.0,0.0,0.0) for v in V2]
        
        plateDim  = [3.0, 3.0]
        pointSpacing = 0.25   
        F3,V3 = triplate(plateDim, pointSpacing)
        V3 = [v.+Point{3,Float64}(1.0+boxDim[1]+plateDim[1]/2.0+2.0,0.0,0.0) for v in V3]

        F, V, C = joingeom([F1, F2, F3], [V1, V2, V3])

        plot_joingeom(F,V,C, testCase)
    end
end