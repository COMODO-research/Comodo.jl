using Comodo
using GLMakie
using GeometryBasics
using FileIO
using Statistics
using Rotations

# Example geometry
testCase = 3
if testCase == 1
    F,V = geosphere(2,1.0)
elseif testCase == 2
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = tofaces(faces(M))
    V = topoints(coordinates(M))        
    F,V,ind1,ind2 = mergevertices(F,V)
    # F,V=subtri(F,V,1)
elseif testCase == 3
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = tofaces(faces(M))
    V = topoints(coordinates(M))    
    F,V,ind1,ind2 = mergevertices(F,V)
end

# Input parameters
p = mean(V,dims=1)[1]; # Point on cutting plane
n = normalizevector(Vec{3, Float64}(0.0,1.0,1.0))# Cutting plane normal
snapTolerance = 1e-6

cutType = :full
Fn,Vn,Cn = trisurfslice(F,V,n,p; output_type=cutType)
Fn,Vn = separate_vertices(Fn,Vn)
CnV = simplex2vertexdata(Fn,Cn)
Mn = GeometryBasics.Mesh(Vn,Fn)


## Visualization
cmap = colormap = cgrad(:Spectral, 5, categorical = true)

s = 1.25*maximum([maximum(map(v-> v[i],V)) - minimum(map(v-> v[i],V)) for i ∈ 1:3])

R = rotation_between(n,[0.0,0.0,1.0])
plateDim = (s,s)
plateElem = (1,1)
FG1,VG1 = quadplate(plateDim,plateElem)
VGn = [GeometryBasics.Point{3, Float64}(R'*v)+p for v ∈ VG1]
MG = GeometryBasics.Mesh(VGn,FG1)

fig = Figure(size=(800,800))

# ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A sliced mesh")
ax1 = LScene(fig[1,1]); 
cc1 = cam3d!(ax1.scene, clipping_mode = :view_relative, projectiontype = Makie.Perspective) 

hp1 = wireframe!(ax1,MG, linewidth=5, color=:red)
hp2 = poly!(ax1,Mn, color=CnV, strokewidth=1, strokecolor=:black, shading = FastShading, transparency=false, colorrange = (-2.5,2.5),colormap=cmap)
hp3 = Colorbar(fig[1,2],hp2,ticks=-2:1:2)

cc1.near[] = 1f-3
cc1.far[] = 100

# ax2 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A sliced mesh")
ax2 = LScene(fig[1,3]); 
cc2 = cam3d!(ax1.scene, clipping_mode = :view_relative, projectiontype = Makie.Perspective) 

Mn = GeometryBasics.Mesh(Vn,Fn[Cn.<=0])
hp4 = poly!(ax2,Mn, color=:white, strokewidth=1, strokecolor=:black, shading = FastShading, transparency=false, colorrange = (-2.5,2.5),colormap=cmap)

cc2.near[] = 1f-3
cc2.far[] = 100



stepRange = range(-s,s,100)
hSlider = Slider(fig[2, :], range = stepRange, startvalue = 0,linewidth=30)

on(hSlider.value) do stepIndex 
    pp = p + stepIndex*n
    Fn,Vn,Cn = trisurfslice(F,V,n,pp; output_type=cutType) 
   
    if isempty(Fn)
        Mn = GeometryBasics.Mesh(V,F)
        CnV = zeros(length(V))
    else
        Fn,Vn = separate_vertices(Fn,Vn)
        CnV = simplex2vertexdata(Fn,Cn)
        Mn = GeometryBasics.Mesh(Vn,Fn)
    end
    
    VGn = [GeometryBasics.Point{3, Float64}(R'*v)+pp for v ∈ VG1] # Rotate plane    
    MG = GeometryBasics.Mesh(VGn,FG1)

    hp1[1] = MG
    hp2[1] = Mn
    hp2.color = CnV

    Mn = GeometryBasics.Mesh(Vn,Fn[Cn.<=0])
    hp4[1] = Mn

end

slidercontrol(hSlider,ax1)
fig