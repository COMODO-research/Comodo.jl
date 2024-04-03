using Comodo
using GLMakie
using GeometryBasics
using FileIO

# Example geometry
testCase = 7
if testCase == 1
    r = 25.25
    F,V = geosphere(5,r)  
    # V = [GeometryBasics.Point{3, Float64}(v[1],v[2],3*v[3]) for v ∈ V]  
elseif testCase==2
    r = 25.25
    F,V = quadsphere(3,r)    
elseif testCase==3
    r = sqrt(3)
    M = cube(r)
    F = faces(M)
    V = coordinates(M)
    # F = quad2tri(F,V; convert_method = "angle")
elseif testCase==4   
    r = 25
    s = r/10
    nc = round(Int64,(2*pi*r)/s)
    d = r*2
    Vc = circlepoints(r, nc; dir=:cw)
    num_steps = round(Int64,d/s)
    num_steps = num_steps + Int64(iseven(num_steps))
    F, V = extrudecurve(Vc, d; s=1, n=[0.0,0.0,1.0], num_steps=num_steps, close_loop=true, face_type=:quad)
elseif testCase==5 # Merged STL for single object
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = togeometrybasics_faces(faces(M))
    V = togeometrybasics_points(coordinates(M))
    F,V,_ = mergevertices(F,V)
    F,V = subtri(F,V,2; method = :loop)
elseif testCase==6 # Merged STL for single object
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = togeometrybasics_faces(faces(M))
    V = togeometrybasics_points(coordinates(M))
    F,V,_ = mergevertices(F,V)
elseif testCase==7
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
s = pointspacingmean(F,V)
cMap = :Spectral
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "1st principal curvature")
hp1 = mesh!(ax1,M, color=K1, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(K1)).*0.1.*(-1,1))
# hp1 = poly!(ax1,M, color=K1, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(K1)).*0.1.*(-1,1),strokecolor=:black,strokewidth=2)
# normalplot(ax1,M)
hpn1 = dirplot(ax1,V,U1; color=:black,linewidth=2,scaleval=s,style=:through)
Colorbar(fig[1, 2],hp1)

ax1 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "2nd principal curvature")
hp2 = mesh!(ax1,M, color=K2, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(K2)).*0.1.*(-1,1))
hpn2 = dirplot(ax1,V,U2; color=:black,linewidth=2,scaleval=s,style=:through)
Colorbar(fig[1, 4],hp1)

ax1 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Mean curvature")
hp2 = mesh!(ax1,M, color=H, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(H)).*0.1.*(-1,1))
Colorbar(fig[2, 2],hp1)

ax1 = Axis3(fig[2, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Gaussian curvature")
hp2 = mesh!(ax1,M, color=G, shading = FastShading, transparency=false,colormap=cMap,colorrange = maximum(abs.(G)).*0.025.*(-1,1))
Colorbar(fig[2, 4],hp1)

fig

