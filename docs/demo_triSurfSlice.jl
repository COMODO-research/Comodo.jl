using Gibbon
using GLMakie
using GeometryBasics
using Statistics
using Rotations
using LinearAlgebra
using FileIO
using SparseArrays

# Example geometry
testCase = 3
if testCase == 1
    F,V = geosphere(2,1.0)
elseif testCase == 2
    # Loading a mesh
    fileName_mesh = joinpath(gibbondir(),"assets","stl","stanford_bunny_low.stl")
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
    fileName_mesh = joinpath(gibbondir(),"assets","stl","david.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = faces(M)
    V = coordinates(M)
    F = togeometrybasics_faces(F)
    V = togeometrybasics_points(V)
    F,V,ind1,ind2 = mergevertices(F,V)
end

# function trisurfslice(F,V; p=(0,0,0))
# Input parameters
p = mean(V,dims=1)[1]; # Point on cutting plane
n = normalizevector(Vec{3, Float64}(0.0,1.0,1.0))# Cutting plane normal
snapTolerance = 1e-6

Fn,Vn = trisurfslice(F,V,n,p; output_type="below")
Mn = GeometryBasics.Mesh(Vn,Fn)


## Visualization
s = 1.25*maximum([maximum(map(v-> v[i],V)) - minimum(map(v-> v[i],V)) for i ∈ 1:3])

R = rotation_between(n,[0.0,0.0,1.0])
plateDim = (s,s)
plateElem = (1,1)
FG1,VG1 = quadplate(plateDim,plateElem)
MG = GeometryBasics.Mesh(VG1,FG1)
fig = Figure(size=(800,800))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A sliced mesh")

stepRange = range(-s/2,s/2,50)
hSlider = Slider(fig[2, 1], range = stepRange, startvalue = 0,linewidth=30)

Mn = lift(hSlider.value) do stepIndex       
    pp = [p[1],p[2],p[3]+stepIndex]
    Fn,Vn = trisurfslice4(F,V,n,pp; output_type="below")        
    return GeometryBasics.Mesh(Vn,Fn)
end

MG = lift(hSlider.value) do stepIndex   
    pp = [p[1],p[2],p[3]+stepIndex]
    
    VGn = [GeometryBasics.Point{3, Float64}(R'*v) for v ∈ VG1] # Rotate plane
    VGn = map(v-> v.+pp,VGn) # Offset plate    
    return GeometryBasics.Mesh(togeometrybasics_points(VGn),FG1)
end

hp1 = mesh!(ax1,GeometryBasics.Mesh(V,F),color=:white, shading = FastShading, transparency=true)
hp2 = wireframe!(ax1,MG, linewidth=5, color=:red)
hp3 = poly!(ax1,Mn, strokewidth=1,color=:white, strokecolor=:blue, shading = FastShading, transparency=false)
# hp3 = normalplot(ax1,Mn)

slidercontrol(hSlider,ax1)
fig
