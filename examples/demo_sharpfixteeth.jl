GLMakie.closeall()
for testCase = 1:2
    if testCase == 1 
        # Create example geometry
        n = 3
        r = 5.0
        F,V = geosphere(n,r)    
        
        # Creating Boolean for faces 
        VF = simplexcenter(F,V)
        B = [vf[3]<r/1.5 for vf in VF]
        
        # "Fix" Boolean
        B_fix1 = sharpfixteeth(F,B; method=:add)
        B_fix2 = sharpfixteeth(F,B; method=:remove)

        # Visualization
        E_cut = boundaryedges(F[B])
        E1 = boundaryedges(F[B_fix1])
        E2 = boundaryedges(F[B_fix2])

        Fp,Vp = separate_vertices(F,V)
        Cp = simplex2vertexdata(Fp,B)
        Cp2 = simplex2vertexdata(Fp,B_fix1)
        Cp3 = simplex2vertexdata(Fp,B_fix2)

        cmap = Makie.Categorical(:viridis) 
        strokewidth = 1
        linewidth = 3

        fig = Figure(size=(1200,800))

        ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
        hp1 = poly!(ax1,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp, shading = FastShading, colormap=cmap)
        hp2 = wireframe!(ax1,GeometryBasics.Mesh(V,E_cut), linewidth=5,color=:red)

        ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
        hp3 = poly!(ax2,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp2, shading = FastShading, colormap=cmap)
        hp4 = wireframe!(ax2,GeometryBasics.Mesh(V,E_cut), linewidth=5,color=:red)
        hp5 = wireframe!(ax2,GeometryBasics.Mesh(V,E1), linewidth=5,color=:cyan)

        ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Removed outward "teeth" """)
        hp6 = poly!(ax3,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp3, shading = FastShading, colormap=cmap)
        hp7 = wireframe!(ax3,GeometryBasics.Mesh(V,E_cut), linewidth=5,color=:red)
        hp8 = wireframe!(ax3,GeometryBasics.Mesh(V,E2), linewidth=5,color=:blue)
        Legend(fig[1, 4], [hp2, hp5, hp8], ["Before","After: added", "After: removed"])
        # Colorbar(fig[1, 5], hp3)
        display(GLMakie.Screen(),fig)

    elseif testCase == 2
        # Create example geometry 
        n = 3
        r = 5.0
        F1,V1 = geosphere(n,r)
        E,V,CE,Fb,Cb = tetgenmesh(F1,V1) 

        # Create Boolean for elements
        VE = simplexcenter(E,V)    
        B = [v[3] < r/1.5 for v in VE]
    
        # "fix" Boolean
        B_fix1 = sharpfixteeth(E,B; method=:add)
        B_fix2 = sharpfixteeth(E,B; method=:remove)

        # Visualisation
        F = element2faces(E) # Triangular faces  
        B_F = repeat(B,inner=4) # Create equivalent for faces
        B_fix1_F = repeat(B_fix1,inner=4) # Create equivalent for faces
        B_fix2_F = repeat(B_fix2,inner=4) # Create equivalent for faces

        indBoundaryFaces = boundaryfaceindices(F)
        Fb = F[indBoundaryFaces]
        B_Fb = B_F[indBoundaryFaces]
        B_fix1_Fb = B_fix1_F[indBoundaryFaces]
        B_fix2_Fb = B_fix2_F[indBoundaryFaces]

        F_cut = element2faces(E[B]) # Triangular faces  
        F_cutb = boundaryfaces(F_cut)

        F1 = element2faces(E[B_fix1])
        F1b = boundaryfaces(F1)

        F2 = element2faces(E[B_fix2])
        F2b = boundaryfaces(F2)

        Fp,Vp = separate_vertices(Fb,V)
        Cp = simplex2vertexdata(Fp,B_Fb)
        Cp2 = simplex2vertexdata(Fp,B_fix1_Fb)
        Cp3 = simplex2vertexdata(Fp,B_fix2_Fb)

        #Visualization
        cmap = cgrad(:viridis, 2, categorical = true)
        strokewidth = 1
        linewidth = 3



        fig = Figure(size=(1200,800))

        ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
        hp1 = poly!(ax1,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp, shading = FastShading, colormap=cmap)
        ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
        hp2 = poly!(ax2,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp2, shading = FastShading, colormap=cmap)
        ax3 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Removed outward "teeth" """)
        hp5 = poly!(ax3,GeometryBasics.Mesh(Vp,Fp), strokewidth = strokewidth, color = Cp3, shading = FastShading, colormap=cmap)
        ax4 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
        hp4 = poly!(ax4,GeometryBasics.Mesh(V,F_cutb), strokewidth = strokewidth, color =:white, shading = FastShading, colormap=cmap)
        ax5 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
        hp6 = poly!(ax5,GeometryBasics.Mesh(V,F1b), strokewidth = strokewidth, color =:white, shading = FastShading, colormap=cmap)
        ax6 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Removed outward "teeth" """)
        hp8 = poly!(ax6,GeometryBasics.Mesh(V,F2b), strokewidth = strokewidth, color =:white, shading = FastShading, colormap=cmap)
        display(GLMakie.Screen(),fig)

    end 
end
