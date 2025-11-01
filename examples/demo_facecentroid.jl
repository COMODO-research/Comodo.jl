using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.Statistics

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1
        r = 1.0
        F,V = geosphere(2,r)   
    elseif testCase == 2    
        r = 1.0    
        F,V = dodecahedron(r)  
    elseif testCase == 3
        r = 1.0
        Fs,Vs = geosphere(2,r) 
        F,V = meshdual(Fs,Vs)
    end
    
    # Compute polygon centroid
    VC = facecentroid(F, V)

    ## Visualization
    fig = Figure(size=(900,900))
    ax1 = AxisGeom(fig[1, 1])
    hp1 = meshplot!(ax1, F, V; color=(:white, 0.9), transparency=true, strokewidth=2.0)
    hp2 = scatter!(ax1, VC, markersize=25, color = :red)
 
    Legend(fig[1, 2], [hp1, hp2],["Faces", "Centroids"])
    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end