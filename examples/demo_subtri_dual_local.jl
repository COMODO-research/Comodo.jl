using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Statistics

#=
This demo shows the use of `subtri_dual_local` which uses the `subtri_dual` 
function to locally refine triangulated meshes. 
=#

GLMakie.closeall()

for testCase = 1:2
    if testCase == 1
        r = 160.0 # Radius
        n1 = 3
        F1, V1 = tridisc(r,n1)

        d = 100.0
        B1 = [norm(mean(V1[f]))<d for f in F1]

        smooth = true        
        F2, V2 = subtri_dual_local(F1, V1, B1; smooth=smooth)

        d = 50.0
        B2 = [norm(mean(V2[f]))<d for f in F2]
        F3, V3= subtri_dual_local(F2, V2, B2)

        ## Visualization
        strokewidth1 = 1
        lineWidth = 4

        fig = Figure(size=(1600,800))

        ax1 = AxisGeom(fig[1, 1], title = "Original", azimuth=-pi/2, elevation=pi/2)
        meshplot!(ax1, F1, V1)

        ax2 = AxisGeom(fig[1, 2], title = "Once refined", azimuth=-pi/2, elevation=pi/2)
        meshplot!(ax2, F2, V2)

        ax3 = AxisGeom(fig[1, 3], title = "Twice refined", azimuth=-pi/2, elevation=pi/2)
        meshplot!(ax3, F3, V3)

        screen = display(GLMakie.Screen(), fig)
    elseif testCase == 2 
        r = 10.0 # Radius
        n1 = 3
        F1, V1 = geosphere(n1, r)

        B1 = [mean(V1[f])[1]< -r/2.0 for f in F1]

        smooth = true        
        F2, V2 = subtri_dual_local(F1, V1, B1; smooth=smooth)

        ## Visualization
        strokewidth1 = 1
        lineWidth = 4

        fig = Figure(size=(1600,800))

        ax1 = AxisGeom(fig[1, 1], title = "Original")
        meshplot!(ax1, F1, V1)

        ax2 = AxisGeom(fig[1, 2], title = "Once refined")
        meshplot!(ax2, F2, V2)

        screen = display(GLMakie.Screen(), fig)
    end
end