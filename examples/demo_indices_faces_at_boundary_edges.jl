using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of `boundaryfaceindices` to obtain the indices of boundary faces of a 
volumetric mesh. 
=#

GLMakie.closeall()

for testCase = 1:5

    if testCase == 1 
        # Square
        w = 1.0
        f,V = square(w)    
        F = quad2tri([f],V)
    elseif testCase == 2 
        # Test quad ring 
        pointSpacing = 0.5
        r1 = 1.0
        r2 = r1-pointSpacing
        n = ceil(Int,((2*pi)*r1)./pointSpacing)
        V1 = circlepoints(r1,10)
        V2 = [v*(r2/r1) for v in V1]        
        F,V = loftlinear(V1,V2; num_steps=2, close_loop=true, face_type=:quad)                        
    elseif testCase == 3
        # Test quad ring refined
        pointSpacing = 0.5
        r1 = 1.0
        r2 = r1-pointSpacing
        n = ceil(Int,((2*pi)*r1)./pointSpacing)
        V1 = circlepoints(r1,10)
        V2 = [v*(r2/r1) for v in V1]        
        F,V = loftlinear(V1,V2; num_steps=2, close_loop=true, face_type=:quad)  
        F,V = subquad(F,V,2; method=:Catmull_Clark)
    elseif testCase == 4            
        # Test tri ring 
        pointSpacing = 0.5
        r1 = 1.0
        r2 = r1-pointSpacing
        n = ceil(Int,((2*pi)*r1)./pointSpacing)
        V1 = circlepoints(r1,10)
        V2 = [v*(r2/r1) for v in V1]        
        F,V = loftlinear(V1,V2; num_steps=2, close_loop=true, face_type=:forwardslash)
    elseif testCase == 5
        # Test tri ring refined
        pointSpacing = 0.5
        r1 = 1.0
        r2 = r1-pointSpacing
        n = ceil(Int,((2*pi)*r1)./pointSpacing)
        V1 = circlepoints(r1,10)
        V2 = [v*(r2/r1) for v in V1]        
        F,V = loftlinear(V1,V2; num_steps=2, close_loop=true, face_type=:forwardslash)        
        F, V = subtri(F,V,2; method=:Loop)
    end

    ind = indices_faces_at_boundary_edges(F)

    # Visualisation
    fig = Figure(size=(1600,800))

    ax1 = AxisGeom(fig[1, 1]; title = "Boundary faces")
    hp2 = edgeplot!(ax1, F, V; color = :black)
    hp3 = meshplot!(ax1, F[ind], V; color = :red)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end