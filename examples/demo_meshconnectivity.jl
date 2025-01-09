using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

testCase = 1
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
E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    

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

ax=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "edge-face")
hp=poly!(ax,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white, shading=FastShading, overdraw=false)

hp1=wireframe!(ax,GeometryBasics.Mesh(V,[E_uni[indE]]),linewidth=5, transparency=true, depth_shift=-1.0f-3, color=:blue)
hp2=poly!(ax,GeometryBasics.Mesh(V,F[con_E2F[indE]]), strokewidth=4,color=:red, shading=FastShading, overdraw=false)

ax=Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face-edge")
hp=poly!(ax,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white, shading=FastShading, overdraw=false)

hp1=poly!(ax,GeometryBasics.Mesh(V,[F[indF]]), strokewidth=4,color=:blue, shading=FastShading, overdraw=false)
hp2=wireframe!(ax,GeometryBasics.Mesh(V,E_uni[con_F2E[indF]]),linewidth=5, transparency=true, depth_shift=-1.0f-3, color=:red)

ax=Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face-face (via edges)")
hp=poly!(ax,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white, shading=FastShading, overdraw=false)

hp1=poly!(ax,GeometryBasics.Mesh(V,[F[indF]]), strokewidth=4,color=:blue, shading=FastShading, overdraw=false)
hp2=poly!(ax,GeometryBasics.Mesh(V,F[con_F2F[indF]]), strokewidth=4,color=:red, shading=FastShading, overdraw=false)

ax=Axis3(fig[1, 4], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "vertex-face")
hp=poly!(ax,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white, shading=FastShading, overdraw=false)

hp1=scatter!(ax,V[indV],color =:blue, markersize = 25)
hp2=poly!(ax,GeometryBasics.Mesh(V,F[con_V2F[indV]]), strokewidth=4,color=:red, shading=FastShading, overdraw=false)

ax=Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "vertex-edge")
hp=poly!(ax,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white, shading=FastShading, overdraw=false)

hp1=scatter!(ax,V[indV],color =:blue, markersize = 25)
hp2=wireframe!(ax,GeometryBasics.Mesh(V,E_uni[con_V2E[indV]]),linewidth=5, transparency=true, depth_shift=-1.0f-3, color=:red)


ax=Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "vertex-vertex")
hp=poly!(ax,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white, shading=FastShading, overdraw=false)

hp1=scatter!(ax,V[indV],color =:blue, markersize = 25)
hp2=scatter!(ax,V[con_V2V[indV]],color =:red, markersize = 25)

ax=Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "edge-edge")
hp=poly!(ax,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white, shading=FastShading, overdraw=false)

hp1=wireframe!(ax,GeometryBasics.Mesh(V,[E_uni[indE]]), linewidth=5, transparency=true, depth_shift=-1.0f-3, color=:blue)
hp2=wireframe!(ax,GeometryBasics.Mesh(V,E_uni[con_E2E[indE]]),linewidth=5, transparency=true, depth_shift=-1.0f-3, color=:red)

ax=Axis3(fig[2, 4], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "face-face (via vertices)")
hp=poly!(ax,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white, shading=FastShading, overdraw=false)

hp1=poly!(ax,GeometryBasics.Mesh(V,[F[indF]]), strokewidth=4,color=:blue, shading=FastShading, overdraw=false)
hp2=poly!(ax,GeometryBasics.Mesh(V,F[con_F2F_v[indF]]), strokewidth=4,color=:red, shading=FastShading, overdraw=false)

fig