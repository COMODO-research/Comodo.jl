using Comodo
using GLMakie
using GeometryBasics
using FileIO
using LinearAlgebra

testCase = 2

if testCase==1 
    # A "zig-zag" curve with known angles for testing
    r = 2.0
    B = [170.0,150.0,135.0,90.0,75.0,60.0,45.0,-45.0,-60.0,-75,-90.0,-135.0,-150.0,-170.0]
    V = Vector{Point{3, Float64}}()
    push!(V,[0.0, 0.0, 0.0])
    for b in B
        push!(V,V[end]+[r, 0.0, 0.0])
        push!(V,V[end]+[r*cosd(b),r*sind(b), 0.0])
    end
    push!(V,V[end]+[r, 0.0, 0.0])

    F,V = extrudecurve(V; extent=1, direction=:positive, n=Vec{3, Float64}(0.0,0.0,-1.0), close_loop=false,face_type=:quad)
elseif testCase==2 
    # An imported STL based geometry of an engineering part with various flat faces and some known angles (e.g. 0, 45, and 90 degrees)
    # This model is relatively evenly sampled and of a relatively low resolution
    fileName_mesh = joinpath(comododir(),"assets","stl","spur_gear_01.stl")
    M = load(fileName_mesh)
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    F,V,_,_ = mergevertices(F,V)
elseif testCase==3 # Same as 2 but higher angular/general resolution
    # An imported STL based geometry of an engineering part with various flat faces and some known angles (e.g. 0, 45, and 90 degrees)
    # This model features a non-homogeneous mesh (includes sharp triangles) and is of a relatively high resolution
    fileName_mesh = joinpath(comododir(),"assets","stl","spur_gear_02.stl")
    M = load(fileName_mesh)
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    F,V,_,_ = mergevertices(F,V)
end

# E = meshedges(F)
# E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)        
# con_E2F = con_edge_face(F,E_uni)
# con_F2E = con_face_edge(F,E_uni,indReverse)   

A,E,con_E2F = edgefaceangles(F,V; deg=true)


G = faceanglesegment(F,V; deg=true, angleThreshold = 22.5, indStart = 1)

## Visualization
linewidth = 3
c = cgrad(:Spectral,maximum(G),categorical = true)
A,E,con_E2F = edgefaceangles(F,V; deg=true)

Lp = .~isnan.(A) .&& abs.(A).>30
Ap = A[Lp]
Ep = E[Lp]

En,Vn_E = separate_vertices(Ep,V)
An = simplex2vertexdata(En,Ap,Vn_E)


Fn,Vn = separate_vertices(F,V)
Gn = simplex2vertexdata(Fn,G,Vn)## Visualization
linewidth = 3

fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Imported mesh")
hp1 = poly!(ax1,GeometryBasics.Mesh(Vn,Fn), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=0.25)
# normalplot(ax1,F,V)
# hp_A = wireframe!(ax1,GeometryBasics.Mesh(Vn_E,En),linewidth=linewidth, transparency=false, color=An,colormap=:Spectral,colorrange = (-120, 120))
# Colorbar(fig[1, 2],hp_A, label = "Angles")

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Groups")
hp2 = poly!(ax2,GeometryBasics.Mesh(Vn,Fn), color=Gn, shading = FastShading, transparency=false,colormap=c)#,strokecolor=:black,strokewidth=0.25)
# normalplot(ax1,F,V)

hp3 = wireframe!(ax2,GeometryBasics.Mesh(Vn_E,En),linewidth=linewidth, transparency=false, color=:black)

Colorbar(fig[1, 3],hp2, label = "Groups")

fig
