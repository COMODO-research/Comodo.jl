using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

GLMakie.closeall()

for testCase = 1:2
    if testCase == 1    
        r = 1.0
        n = 1
        Fb,Vb = geosphere(n,r)
        E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb)  
    elseif testCase == 2
        nSmooth = 50
        r = 5.0
        f = 0.75
        pointSpacing = 1.0
        E,V = hexsphere(r,pointSpacing; f=f, nSmooth=nSmooth) 
    end

    edgeSet = elementEdges(E)

    F = element2faces(E) # Triangular faces

    ## Visualization
    fig = Figure(size=(800,800))

    ax1 = AxisGeom(fig[1, 1], title = "Mesh edges")
    hp1 = meshplot!(ax1, F, V, strokewidth=0.0)
    hp2 = edgeplot!(ax1, edgeSet, V, color=:red, linewidth=3)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end