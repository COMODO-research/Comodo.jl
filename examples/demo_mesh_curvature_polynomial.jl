using Comodo
using GLMakie
using GeometryBasics
using FileIO
using LinearAlgebra
using Statistics
using Rotations

# Example geometry
testCase = 6
if testCase == 1
    r = 25.25
    F,V = geosphere(5,r)  
    # V = [GeometryBasics.Point{3, Float64}(v[1],v[2],3*v[3]) for v ∈ V]  
elseif testCase==2
    r = 1.0
    F,V = quadsphere(3,r)    
elseif testCase==3
    r=sqrt(3)
    M = cube(r)
    F=faces(M)
    V=coordinates(M)
    # F = quad2tri(F,V; convert_method = "angle")
elseif testCase==4 # Merged STL for single object
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = togeometrybasics_faces(faces(M))
    V = togeometrybasics_points(coordinates(M))
    F,V,_ = mergevertices(F,V)
    F,V = subtri(F,V,2; method = :loop)
elseif testCase==5 # Merged STL for single object
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = togeometrybasics_faces(faces(M))
    V = togeometrybasics_points(coordinates(M))
    F,V,_ = mergevertices(F,V)
elseif testCase==6
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","obj","motherChild_5k.obj")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = togeometrybasics_faces(faces(M))
    V = togeometrybasics_points(coordinates(M))

end
    
M = GeometryBasics.Mesh(V,F)

K1,K2,U1,U2,H,G = mesh_curvature_polynomial(F,V)

## Visualization

cMap = :Spectral
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "1st principal curvature")
hp1 =mesh!(ax1,M, color=K1, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(K1)).*0.1.*(-1,1))
hpn1 = dirplot(ax1,V,U1; color=:black,linewidth=1,scaleval=3,style=:through)
Colorbar(fig[1, 2],hp1)

ax1 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "2nd principal curvature")
hp2 =mesh!(ax1,M, color=K2, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(K2)).*0.1.*(-1,1))
hpn2 = dirplot(ax1,V,U2; color=:black,linewidth=1,scaleval=3,style=:through)
Colorbar(fig[1, 4],hp1)

ax1 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Mean curvature")
hp2 =mesh!(ax1,M, color=H, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(H)).*0.1.*(-1,1))
Colorbar(fig[2, 2],hp1)

ax1 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Gaussian curvature")
hp2 =mesh!(ax1,M, color=G, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(G)).*0.025.*(-1,1))
Colorbar(fig[2, 4],hp1)

fig


