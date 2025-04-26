using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of `boundaryfaces` to obtain the boundary faces of a 
volumetric mesh. 
=#

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1 
        sampleSize = 10
        pointSpacing = 2.0
        boxDim = sampleSize.*[1.0,1.0,1.0] # Dimensionsions for the box in each direction
        boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 
        E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)    

        # Get boundary faces from either the total set of element faces or the elements
        # Fb = boundaryfaces(F)
        Fb = boundaryfaces(E)
        Cb = ones(length(F))
    elseif testCase == 2 
        sampleSize = 10
        pointSpacing = 2.0
        boxDim = sampleSize.*[1.0,1.0,1.0] # Dimensionsions for the box in each direction
        boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 
        E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl) 
        VE = simplexcenter(E,V)
        elementLabels = [v[3]<=eps(0.0) for v in VE]

        # Get boundary faces from either the total set of element faces or the elements    
        Fb = boundaryfaces(E; elementLabels=elementLabels)
        Cb = ones(length(F))
    elseif testCase == 3
        F1,V1 = geosphere(3,1.0)
        E,V,CE,Fb,Cb = tetgenmesh(F1,V1)   
        F = element2faces(E)

        # Get boundary faces from either the total set of element faces or the elements
        Fb = boundaryfaces(F)    
        Cb = ones(length(F))
    end

    # Visualisation
    fig = Figure(size=(1600,800))

    Fs,Vs = separate_vertices(F,V)
    Fbs,Vbs = separate_vertices(Fb,V)
    Cbs = simplex2vertexdata(Fbs,Cb)

    ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Input faces")
    hp1 = poly!(ax1,GeometryBasics.Mesh(Vs,Fs), strokewidth=1,shading=FastShading,strokecolor=:black, color=:white, transparency=true, overdraw=false)

    ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary faces")
    hp1 = poly!(ax2,GeometryBasics.Mesh(Vbs,Fbs), strokewidth=1,shading=FastShading,strokecolor=:black, color=Cbs, transparency=true, overdraw=false)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end