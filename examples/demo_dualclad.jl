using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.GLMakie.Colors
using FileIO

c1 = RGBf(1.0, 0.30196078431372547, 0.023529411764705882)
c2 = RGBf(0.2235294117647059, 1.0, 0.0784313725490196)
c3 = RGBf(0.5803921568627451, 0.3411764705882353, 0.9215686274509803)

#=
This demo shows the use of the dualclag function. 
=#

testCase = 1
if testCase == 1 
    F,V = geosphere(2,1.0)
elseif testCase == 2 
    F,V = geosphere(1,1.0)
    B = [v[3]>0 for v in V]
    BF = [all(B[f]) for f in F]
    F = F[BF]
    F,V = remove_unused_vertices(F,V)
elseif testCase == 3
    F,V = tetrahedron(√3)        
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
con_type = :face
Fs,Fq,Vs = dualclad(F,V,s; connectivity=con_type)

# Visualisation
strokewidth =1

Vn = V-0.01*vertexnormal(F,V)

fig = Figure(size = (1200,1200))
ax = Axis3(fig[1, 1], aspect = :data)

hp1 = mesh!(ax, GeometryBasics.Mesh(Vn,F), color=c3,transparency=false, shading = FastShading)
hp2 = poly!(ax, GeometryBasics.Mesh(Vs,Fs), color=c1,transparency=false,strokewidth=strokewidth,strokecolor=:white, shading = FastShading)
hp3 = poly!(ax, GeometryBasics.Mesh(Vs,Fq), color=c2,transparency=false,strokewidth=strokewidth,strokecolor=:white, shading = FastShading)

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

# fileName = comododir()*"/assets/temp/dualclad_anim_06.mp4"
# slider2anim(fig,hSlider,fileName; backforth=true, duration=6)

fig


