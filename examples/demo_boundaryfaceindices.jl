using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.GLMakie.Colors

#=
This demo shows the use of `boundaryfaceindices` to obtain the indices of boundary faces of a 
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

        # Get boundary face indices 
        indBoundaryFaces = boundaryfaceindices(F)

        # Use indices to obtain boundary faces
        Fb = F[indBoundaryFaces]
        Cb = ones(length(Fb))
    elseif testCase == 2 
        sampleSize = 10
        pointSpacing = 2.0
        boxDim = sampleSize.*[1.0,1.0,1.0] # Dimensionsions for the box in each direction
        boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 
        E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)    
        VE = simplexcenter(E,V)
        elementLabels = [v[3]<=eps(0.0) for v in VE]

        # Get boundary face indices 
        indBoundaryFaces = boundaryfaceindices(F; elementLabels=elementLabels)

        # Use indices to obtain boundary faces
        Fb = F[indBoundaryFaces]
        Cb = ones(length(Fb))
    elseif testCase == 3
        F1,V1 = geosphere(3,1.0)
        E,V,CE,Fb,Cb = tetgenmesh(F1,V1)   
        F = element2faces(E)

        # Get boundary face indices 
        indBoundaryFaces = boundaryfaceindices(F)

        # Use indices to obtain boundary faces
        Fb = F[indBoundaryFaces]
        Cb = ones(length(Fb))
    end

    # Visualisation
    fig = Figure(size=(1600,800))

    Fs,Vs = separate_vertices(F,V)
    Fbs,Vbs = separate_vertices(Fb,V)
    Cbs = simplex2vertexdata(Fbs,Cb)

    ax1 = AxisGeom(fig[1, 1]; title = "Input faces")
    hp1 = meshplot!(ax1, Fs, Vs; color = (:white, 0.25), transparency=true)

    ax2 = AxisGeom(fig[1, 2]; title = "Boundary faces")
    hp2 = meshplot!(ax2, Fbs, Vbs; color = Cbs, transparency=true)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end