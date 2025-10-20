using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.GLMakie.Colors
using Comodo.Statistics

#=
This demo shows the use of `hex2tet` to convert hexahedral elements to 
tetrahedral elements. 
=#

GLMakie.closeall()
for testCase = 1:3
    if testCase == 1
        # Example hexahedral mesh
        e, V = hexahedronElement(2.0)
        E = [e]        
        F = element2faces(E)

        # Convert to tetrahedra
        E_tet = hex2tet(E, 1)    

        # Visualization
        fig = Figure(size=(1700, 1000))
        Fs, Vs = separate_vertices(F, V)
        E_tet_s, V_tet_s = separate_vertices(E_tet, V)
        F_tet_s = element2faces(E_tet_s)
        F_tet_s, V_tet_s = separate_vertices(F_tet_s, V_tet_s)

        ax1 = AxisGeom(fig[1,1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", 
        title = "Type 1")
        hp1 = meshplot!(ax1, F_tet_s, V_tet_s, color=(:white, 0.25), strokewidth=3, strokecolor=:blue, transparency=true)
        # hp3 = normalplot(ax1, F_tet_s, V_tet_s)

        E_tet_s, V_tet_s = separate_vertices(E_tet, V; scaleFactor=0.75)
        F_tet_s = element2faces(E_tet_s)
        F_tet_s, V_tet_s = separate_vertices(F_tet_s, V_tet_s)

        ax1 = AxisGeom(fig[1,2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", 
        title = "Type 1")
        hp2 = meshplot!(ax1, F_tet_s, V_tet_s, color=(:red, 1.0), strokewidth=3)
        # hp3 = normalplot(ax1, F_tet_s, V_tet_s)

        screen = display(GLMakie.Screen(), fig)
        GLMakie.set_title!(screen, "testCase = $testCase")
    elseif testCase == 2
        # Example hexahedral mesh
        e, V = hexahedronElement(2.0)
        E = [e]        
        F = element2faces(E)

        # Visualization
        fig = Figure(size=(1700, 1000))
        Fs, Vs = separate_vertices(F, V)

        plotSet = 1:14 
        n1 = ceil(Int,sqrt(length(plotSet))/1.5)
        n2 = ceil(Int, length(plotSet)/n1)
        cartInd = CartesianIndices((n1,n2))
        for (i_plot, i) in enumerate(plotSet)
            E_tet = hex2tet(E, i)    # Get tetrahedral elements for the hex

            E_tet_s, V_tet_s = separate_vertices(E_tet, V; scaleFactor=0.65)
            F_tet_s = element2faces(E_tet_s)
            F_tet_s, V_tet_s = separate_vertices(F_tet_s, V_tet_s)

            ij = cartInd[i_plot]
            ax1 = AxisGeom(fig[ij[1], ij[2]], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", 
            title = "Type $i")
            # hp1 = meshplot!(ax1, Fs, Vs, color=(:white, 0.2), transparency = true)
            if i<3
                c = :blue
            else
                c = :red
            end
            hp2 = meshplot!(ax1, F_tet_s, V_tet_s, color=(c, 1.0), strokewidth=3)
            # hp3 = normalplot(ax1, F_tet_s, V_tet_s)
        end
        screen = display(GLMakie.Screen(), fig)
        GLMakie.set_title!(screen, "testCase = $testCase")
    
    elseif testCase == 3
        # Example hexahedral mesh
        pointSpacing = 0.5
        boxDim = [1.0,1.0,1.0] # Dimensionsions for the box in each direction
        boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 
        E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)

        # Convert to tetrahedra
        meshType = [1, 1, 2, 2, 3, 4, 5, 14]
        E_tet = hex2tet(E, meshType)    

        # Visualization
        fig = Figure(size=(1700, 1000))
        Fs, Vs = separate_vertices(F, V)
        E_tet_s, V_tet_s = separate_vertices(E_tet, V)
        F_tet_s = element2faces(E_tet_s)
        F_tet_s, V_tet_s = separate_vertices(F_tet_s, V_tet_s)

        ax1 = AxisGeom(fig[1,1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", 
        title = "Mixed types")
        hp1 = meshplot!(ax1, F_tet_s, V_tet_s, color=(:white, 0.25), strokewidth=3, strokecolor=:blue, transparency=true)
        # hp3 = normalplot(ax1, F_tet_s, V_tet_s)

        E_tet_s, V_tet_s = separate_vertices(E_tet, V; scaleFactor=0.75)
        F_tet_s = element2faces(E_tet_s)
        F_tet_s, V_tet_s = separate_vertices(F_tet_s, V_tet_s)

        ax1 = AxisGeom(fig[1,2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", 
        title = "Mixed types")
        hp2 = meshplot!(ax1, F_tet_s, V_tet_s, color=(:red, 1.0), strokewidth=3)
        # hp3 = normalplot(ax1, F_tet_s, V_tet_s)


        screen = display(GLMakie.Screen(), fig)
        GLMakie.set_title!(screen, "testCase = $testCase")
    end    
end
