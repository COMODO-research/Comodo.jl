using Comodo
using GLMakie
using GeometryBasics
using FileIO
using StaticArrays

#=
This demo shows the use of the dualclag function. 
=#

testCase = 2
if testCase == 1 
    F,V = geosphere(2,1.0)
elseif testCase == 2 
    F,V = geosphere(1,1.0)
    B = [v[3]>0 for v in V]
    BF = [all(B[f]) for f in F]
    F = F[BF]
    F,V = remove_unused_vertices(F,V)
elseif testCase == 3
    M = tetrahedron(√3)
    F = faces(M)
    V = coordinates(M)
    # n = 2
    # F,V = subquad(F,V,n; method=:Catmull_Clark)
elseif testCase == 4
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    F,V,_ = mergevertices(F,V)
elseif testCase == 5
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    F,V = mergevertices(F,V)
    # n = 1
    # F,V = subtri(F,V,n; method=:Loop)
elseif testCase == 6
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    F,V = mergevertices(F,V)
    B = [v[3]>-10 for v in V]
    BF = [all(B[f]) for f in F]
    F = F[BF]
    F,V = remove_unused_vertices(F,V)
    
end




s = 0.5
con_type = :edge
Fs,Fq,Vs = dualclad(F,V,s; connectivity=con_type)

# Rotate the coordinates
fig = Figure(size = (800,800))
ax = Axis3(fig[1, 1], aspect = :data)

# hp1 = poly!(ax, GeometryBasics.Mesh(V,F), color=:blue,transparency=true,strokewidth=0.15,strokecolor=:black, shading = FastShading)

hp2 = poly!(ax, GeometryBasics.Mesh(Vs,Fs), color=:white,transparency=false,strokewidth=0.15,strokecolor=:black, shading = FastShading)
hp3 = poly!(ax, GeometryBasics.Mesh(Vs,Fq), color=:white,transparency=false,strokewidth=0.15,strokecolor=:black, shading = FastShading)

# hp2 = mesh!(ax, GeometryBasics.Mesh(Vs,Fs), color=:white,transparency=false, shading = FastShading)
# hp3 = mesh!(ax, GeometryBasics.Mesh(Vs,Fq), color=:white,transparency=false, shading = FastShading)

stepRange = range(1.0,0.0,50)
hSlider = Slider(fig[2, :], range = stepRange, startvalue = 1.0,linewidth=30)

slidercontrol(hSlider,fig)

on(hSlider.value) do s
    Fs,Fq,Vs = dualclad(F,V,s; connectivity=con_type)
    
    hp2[1] = GeometryBasics.Mesh(Vs,Fs)
    hp3[1] = GeometryBasics.Mesh(Vs,Fq)
end

set_close_to!(hSlider,0.5)

# fileName = comododir()*"/assets/temp/dualclad_anim_05.mp4"
# slider2anim(fig,hSlider,fileName; backforth=true, duration=3)

fig


