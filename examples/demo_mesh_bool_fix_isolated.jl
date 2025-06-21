using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.LinearAlgebra
using Comodo.Rotations


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
        B_fix1 = deepcopy(B)
        mesh_bool_fix_isolated!(F,B_fix1; method=:add)
        B_fix2 = deepcopy(B)
        mesh_bool_fix_isolated!(F,B_fix2; method=:remove)

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

        ax1 = AxisGeom(fig[1, 1], xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
        hp1 = meshplot!(ax1,Fp, Vp, strokewidth = strokewidth, color = Cp, colormap=cmap)
        hp2 = edgeplot!(ax1,GeometryBasics.Mesh(V,E_cut), linewidth=5,color=:red)

        ax2 = AxisGeom(fig[1, 2], xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
        hp3 = meshplot!(ax2,Fp, Vp, strokewidth = strokewidth, color = Cp2, colormap=cmap)
        hp4 = edgeplot!(ax2,GeometryBasics.Mesh(V,E_cut), linewidth=5,color=:red)
        hp5 = edgeplot!(ax2,GeometryBasics.Mesh(V,E1), linewidth=5,color=:cyan)

        ax3 = AxisGeom(fig[1, 3], xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Removed outward "teeth" """)
        hp6 = meshplot!(ax3,Fp, Vp, strokewidth = strokewidth, color = Cp3, colormap=cmap)
        hp7 = edgeplot!(ax3,GeometryBasics.Mesh(V,E_cut), linewidth=5,color=:red)
        hp8 = edgeplot!(ax3,GeometryBasics.Mesh(V,E2), linewidth=5,color=:blue)
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
        B_fix1 = deepcopy(B)
        mesh_bool_fix_isolated!(E,B_fix1; method=:add)
        B_fix2 = deepcopy(B)
        mesh_bool_fix_isolated!(E,B_fix2; method=:remove)

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

        ax1 = AxisGeom(fig[1, 1],  xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
        hp1 = meshplot!(ax1,Fp, Vp,  strokewidth = strokewidth, color = Cp, colormap=cmap)
        ax2 = AxisGeom(fig[1, 2],  xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
        hp2 = meshplot!(ax2,Fp, Vp,  strokewidth = strokewidth, color = Cp2, colormap=cmap)
        ax3 = AxisGeom(fig[1, 3],  xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Removed outward "teeth" """)
        hp5 = meshplot!(ax3,Fp, Vp,  strokewidth = strokewidth, color = Cp3, colormap=cmap)
        ax4 = AxisGeom(fig[2, 1],  xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original Boolian vector on surface")
        hp4 = meshplot!(ax4,GeometryBasics.Mesh(V,F_cutb), strokewidth = strokewidth, color =:white, colormap=cmap)
        ax5 = AxisGeom(fig[2, 2],  xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Added inward "teeth" """)
        hp6 = meshplot!(ax5,GeometryBasics.Mesh(V,F1b), strokewidth = strokewidth, color =:white, colormap=cmap)
        ax6 = AxisGeom(fig[2, 3],  xlabel = "X", ylabel = "Y", zlabel = "Z", title = """Removed outward "teeth" """)
        hp8 = meshplot!(ax6,GeometryBasics.Mesh(V,F2b), strokewidth = strokewidth, color =:white, colormap=cmap)
        display(GLMakie.Screen(),fig)

    end 
end
