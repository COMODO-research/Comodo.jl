using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

GLMakie.closeall()

for testCase = 1:4
    if testCase == 1
        # Example curves
        nx = 10
        num_steps = 11
        V = gridpoints(range(0.0,5.0,nx),range(0.0,5.0,num_steps),0.0)

        # Periodicity types, first entry is for direction 1 (which will here close over the ends), the second for direction 2 (which will close in the circle direction)
        periodicity1 = (false, false)
        periodicity2 = (false, false)
        periodicity3 = (false, false)
        periodicity4 = (false, false)
        periodicity5 = (false, false)
        periodicity6 = (false, false)   
    elseif testCase == 2
        # Example curves
        r = 1.0
        nc = 2
        t = range(0,2.0*π-(2.0*π/nc),nc)
        Vc = [Point3{Float64}(cos(tt),sin(tt),0.0) for tt ∈ t]
        num_steps = 2
        zLevels = range(0.0,3.0,num_steps)
        V = deepcopy(Vc)
        for i ∈ 2:num_steps
            append!(V, deepcopy(Vc) .+ Point3{Float64}(0.0,0.0,zLevels[i]))
        end
            
        # Periodicity types, first entry is for direction 1 (which will here close over the ends), the second for direction 2 (which will close in the circle direction)
        periodicity1 = (false, false)
        periodicity2 = (false, false)
        periodicity3 = (false, false)
        periodicity4 = (false, false)
        periodicity5 = (false, false)
        periodicity6 = (false, false)
    elseif testCase == 3
        # Example curves
        r = 1.0
        nc = 16
        t = range(0,2.0*π-(2.0*π/nc),nc)
        Vc = [Point3{Float64}(cos(tt),sin(tt),0.0) for tt ∈ t]
        num_steps = 9
        zLevels = range(0.0,3.0,num_steps)
        V = deepcopy(Vc)
        for i ∈ 2:num_steps
            append!(V, deepcopy(Vc) .+ Point3{Float64}(0.0,0.0,zLevels[i]))
        end
            
        # Periodicity types, first entry is for direction 1 (which will here close over the ends), the second for direction 2 (which will close in the circle direction)
        periodicity1 = (false, false)
        periodicity2 = (false, false)
        periodicity3 = (true, false)
        periodicity4 = (true, false)
        periodicity5 = (true, false)
        periodicity6 = (true, false)
    elseif testCase == 4
        # Example curves
        r = 1.0
        nc = 16
        t = range(2.0*π-(2.0*π/nc),0,nc)
        Vc = [Point3{Float64}(2.0+cos(tt),0.0,sin(tt)) for tt ∈ t]
        n = Vec{3, Float64}(0.0,0.0,1.0)
        num_steps = 25

        # Create test data consisting of a revolved circle (resulting in a doughnut shape) 
        _,V = revolvecurve(Vc; extent=(2*pi-(2*pi/num_steps)), direction=:negative, n=n, num_steps=num_steps, periodicity = (false,false),face_type=:quad)

        # Periodicity types, first entry is for direction 1 (which will here close over the ends), the second for direction 2 (which will close in the circle direction)
        periodicity1 = (false, false)
        periodicity2 = (true, false)
        periodicity3 = (false, true)
        periodicity4 = (true, true)
        periodicity5 = (true, true)
        periodicity6 = (true, true)
    end

    VG = deepcopy(V) # Copy original grid as gridpoints manipulates the coordinates for some face types

    # Face types
    face_type1 = :quad
    face_type2 = :forwardslash
    face_type3 = :backslash
    face_type4 = :tri
    face_type5 = :tri_even
    face_type6 = :tri_even

    # Tri method direction 
    tri_dir4 = 1
    tri_dir5 = 1
    tri_dir6 = 2

    F1 = grid2surf(V,num_steps; face_type=face_type1, periodicity=periodicity1)
    F2 = grid2surf(V,num_steps; face_type=face_type2, periodicity=periodicity2)
    F3 = grid2surf(V,num_steps; face_type=face_type3, periodicity=periodicity3)

    V4 = deepcopy(V) # copy since face_type = :tri will modify input vertices
    F4 = grid2surf(V4,num_steps; face_type=face_type4, periodicity=periodicity4, tri_dir=tri_dir4)

    V5 = deepcopy(V) # copy since face_type = :tri will modify input vertices
    F5 = grid2surf(V5,num_steps; face_type=face_type5, periodicity=periodicity5, tri_dir=tri_dir5)

    V6 = deepcopy(V) # copy since face_type = :tri will modify input vertices
    F6 = grid2surf(V6,num_steps; face_type=face_type6, periodicity=periodicity6, tri_dir=tri_dir6)

    E1 = boundaryedges(F1)
    E2 = boundaryedges(F2)
    E3 = boundaryedges(F3)
    E4 = boundaryedges(F4)
    E5 = boundaryedges(F5)
    E6 = boundaryedges(F6)

    M1 = GeometryBasics.Mesh(V,F1)
    M2 = GeometryBasics.Mesh(V,F2)
    M3 = GeometryBasics.Mesh(V,F3)
    M4 = GeometryBasics.Mesh(V4,F4)
    M5 = GeometryBasics.Mesh(V5,F5)
    M6 = GeometryBasics.Mesh(V6,F6)

    ## Visualization
    markerSize = 10
    strokeWidth = 2
    linewidth = 4

    fig = Figure(size=(1600,1200))

    ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """periodicity=$periodicity1, face_type=$face_type1""")
    poly!(ax1,M1, strokewidth=strokeWidth,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
    wireframe!(ax1,GeometryBasics.Mesh(V,E1),linewidth=linewidth, transparency=false, color=:red)
    # normalplot(ax1,F1,V)
    scatter!(ax1,VG,color=:green,markersize=markerSize)

    ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """periodicity=$periodicity2, face_type=$face_type2""")
    poly!(ax2,M2, strokewidth=strokeWidth,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
    wireframe!(ax2,GeometryBasics.Mesh(V,E2),linewidth=linewidth, transparency=false, color=:red)
    scatter!(ax2,VG,color=:green,markersize=markerSize)

    ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """periodicity=$periodicity3, face_type=$face_type3""")
    poly!(ax3,M3, strokewidth=strokeWidth,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
    wireframe!(ax3,GeometryBasics.Mesh(V,E3),linewidth=linewidth, transparency=false, color=:red)
    scatter!(ax3,VG,color=:green,markersize=markerSize)

    ax4 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """periodicity=$periodicity4, face_type=$face_type4, tri_dir=$tri_dir4""")
    poly!(ax4,M4, strokewidth=strokeWidth,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
    wireframe!(ax4,GeometryBasics.Mesh(V4,E4),linewidth=linewidth, transparency=false, color=:red)
    scatter!(ax4,VG,color=:green,markersize=markerSize)

    ax5 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """periodicity=$periodicity5, face_type=$face_type5, tri_dir=$tri_dir5""")
    poly!(ax5,M5, strokewidth=strokeWidth,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
    wireframe!(ax5,GeometryBasics.Mesh(V5,E5),linewidth=linewidth, transparency=false, color=:red)
    scatter!(ax5,VG,color=:green,markersize=markerSize)


    ax6 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """periodicity=$periodicity6, face_type=$face_type6, tri_dir=$tri_dir6""")
    hp2 = poly!(ax6,M6, strokewidth=strokeWidth,color=:white, strokecolor=:black, shading = FastShading, transparency=false)
    hp3 = wireframe!(ax6,GeometryBasics.Mesh(V6,E6),linewidth=linewidth, transparency=false, color=:red)
    hp1 = scatter!(ax6,VG,color=:green,markersize=markerSize)

    Legend(fig[:,4],[hp1,hp2,hp3],["Original input grid", "Output surface", "Boundary"])

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end
