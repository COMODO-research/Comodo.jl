using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

GLMakie.closeall()

for testCase = 1:2
    if testCase==1
        r = 1 
        F,V = platonicsolid(4,r) # Get icosahedron    
        F = F[2:end-1]
        F,V = subtri(F,V,1)
    elseif testCase ==2
        r = 1 
        F,V = platonicsolid(2,r) # Get icosahedron    
        F,V = subquad(F,V,2)
    end

    E = meshedges(F)
    E_uni,indReverse = gunique(E; return_unique=Val(true), return_inverse=Val(true), sort_entries=true)   

    C = meshconnectivity(F,V)

    E_uni = C.edge_vertex

    con_E2E = C.edge_edge
    con_E2F = C.edge_face

    con_F2E = C.face_edge
    con_F2V = C.face_vertex
    con_F2F = C.face_face
    con_F2F_v = C.face_face_v

    con_V2E = C.vertex_edge
    con_V2F = C.vertex_face
    con_V2V = C.vertex_vertex

    indE = 41
    indF = 5
    indV = 2

    fig = Figure(size=(1200,800))

    ax = AxisGeom(fig[1, 1], title = "edge-face")
    hp = meshplot!(ax, F, V, strokewidth=3)

    hp1 = edgeplot!(ax, [E_uni[indE]], V, linewidth=5, color=:blue)
    hp2 = meshplot!(ax, F[con_E2F[indE]], V, strokewidth=4,color=:red)

    ax = AxisGeom(fig[1, 2], title = "face-edge")
    hp = meshplot!(ax, F, V, strokewidth=3)

    hp1 = meshplot!(ax, [F[indF]], V, strokewidth=4,color=:blue)
    hp2 = edgeplot!(ax, E_uni[con_F2E[indF]], V, linewidth=5, color=:red)

    ax = AxisGeom(fig[1, 3], title = "face-face (via edges)")
    hp = meshplot!(ax, F, V, strokewidth=3)

    hp1 = meshplot!(ax, [F[indF]], V, strokewidth=4,color=:blue)
    hp2 = meshplot!(ax, F[con_F2F[indF]], V, strokewidth=4, color=:red)

    ax = AxisGeom(fig[1, 4], title = "vertex-face")
    hp = meshplot!(ax, F, V, strokewidth=3)

    hp1 = scatter!(ax,V[indV],color =:blue, markersize = 25)
    hp2 = meshplot!(ax, F[con_V2F[indV]], V, strokewidth=4, color=:red)

    ax = AxisGeom(fig[2, 1], title = "vertex-edge")
    hp = meshplot!(ax, F, V, strokewidth=3)

    hp1 = scatter!(ax,V[indV],color =:blue, markersize = 25)
    hp2 = edgeplot!(ax, E_uni[con_V2E[indV]], V, linewidth=5, color=:red)


    ax = AxisGeom(fig[2, 2], title = "vertex-vertex")
    hp = meshplot!(ax, F, V, strokewidth=3)

    hp1 = scatter!(ax, V[indV], color =:blue, markersize = 25)
    hp2 = scatter!(ax, V[con_V2V[indV]], color =:red, markersize = 25)

    ax = AxisGeom(fig[2, 3], title = "edge-edge")
    hp = meshplot!(ax, F, V, strokewidth=3)

    hp1 = edgeplot!(ax, [E_uni[indE]], V, linewidth=5,  color=:blue)
    hp2 = edgeplot!(ax, E_uni[con_E2E[indE]], V,linewidth=5,  color=:red)

    ax = AxisGeom(fig[2, 4], title = "face-face (via vertices)")
    hp = meshplot!(ax, F, V, strokewidth=3)

    hp1 = meshplot!(ax, [F[indF]], V, strokewidth=4,color=:blue)
    hp2 = meshplot!(ax, F[con_F2F_v[indF]], V, strokewidth=4,color=:red)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end