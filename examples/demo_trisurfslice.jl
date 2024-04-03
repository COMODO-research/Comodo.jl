using Comodo
using GLMakie
using GeometryBasics
using FileIO
using Statistics
using Rotations

# Example geometry
testCase = 1
if testCase == 1
    F,V = geosphere(2,1.0)
elseif testCase == 2
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = faces(M)
    V = coordinates(M)
    F = togeometrybasics_faces(F)
    V = togeometrybasics_points(V)
    F,V,ind1,ind2 = mergevertices(F,V)
    # F,V=subtri(F,V,1)
elseif testCase == 3
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = faces(M)
    V = coordinates(M)
    F = togeometrybasics_faces(F)
    V = togeometrybasics_points(V)
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
s = 1.25*maximum([maximum(map(v-> v[i],V)) - minimum(map(v-> v[i],V)) for i ∈ 1:3])

R = rotation_between(n,[0.0,0.0,1.0])
plateDim = (s,s)
plateElem = (1,1)
FG1,VG1 = quadplate(plateDim,plateElem)
VGn = [GeometryBasics.Point{3, Float64}(R'*v)+p for v ∈ VG1]
MG = GeometryBasics.Mesh(VGn,FG1)

fig = Figure(size=(800,800))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A sliced mesh")

stepRange = range(-s,s,50)
hSlider = Slider(fig[2, 1], range = stepRange, startvalue = 0,linewidth=30)

# hp1 = mesh!(ax1,GeometryBasics.Mesh(V,F),color=:white, shading = FastShading, transparency=true)
hp2 = wireframe!(ax1,MG, linewidth=5, color=:red)
hp3 = poly!(ax1,Mn, color=CnV, strokewidth=1, strokecolor=:black, shading = FastShading, transparency=false, colorrange = (-2,2),colormap=:Spectral)
# hp3 = normalplot(ax1,Mn)
hp4 = Colorbar(fig[1,2],hp3)

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
    MG = GeometryBasics.Mesh(togeometrybasics_points(VGn),FG1)

    hp2[1] = MG
    hp3[1] = Mn
    hp3.color = CnV
end

slidercontrol(hSlider,ax1)
fig