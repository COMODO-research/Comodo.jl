using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

GLMakie.closeall()

for testCase = 1:2
    if testCase == 1 
        plateDim1 = [10.0,14.0]
        pointSpacing1 = 2.0

        orientation1 = :up
        F1,V1 = triplate(plateDim1,pointSpacing1; orientation=orientation1)
        t = 6
        n = 4
        direction=:positive
        E, V = extrudefaces(F1,V1; extent=t, direction=direction, num_steps=n)
        F = element2faces(E) # Triangular faces
    elseif testCase == 2 
        F1,V1 = hemisphere(2,1.0; face_type=:tri, closed=false)
        t = 0.5
        n = 5
        direction=:negative
        E, V = extrudefaces(F1,V1; extent=t, direction=direction, num_steps=n)
        F = element2faces(E) # Triangular faces
    end

    E_penta15, V_penta15 = penta6_penta15(E,V)
    F_penta15 = element2faces(E_penta15)

    ## Visualization
    cmap = cgrad(:Spectral, 5, categorical = true)

    strokewidth = 1 

    fig = Figure(size=(800,800))

    ax1 = AxisGeom(fig[1, 1], title = "Penta5 faces")
    hp1 = meshplot!(ax1, F[1], V)
    hp2 = meshplot!(ax1, F[2], V)
    scatter!(ax1,V,color=:black,markersize=10)

    ax2 = AxisGeom(fig[1, 2], title = "Penta15 faces")
    hp3 = meshplot!(ax2, F_penta15[1], V_penta15)
    hp4 = meshplot!(ax2, F_penta15[2], V_penta15)
    scatter!(ax2,V_penta15,color=:black,markersize=10)
    # normalplot(ax2,F_penta15[2],V_penta15)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end