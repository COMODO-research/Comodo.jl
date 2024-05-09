using Comodo
using GLMakie
using GeometryBasics
using FileIO
using Rotations

#=
This demo shows the use of the `scalesimplex` function. 
=#

# Loading a mesh

# fileName_mesh = joinpath(comododir(),"assets","obj","spot_control_mesh_texture.obj")
# M = load(fileName_mesh)

fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
M = load(fileName_mesh)

# Obtain mesh faces and vertices
F = tofaces(faces(M))
V = topoints(coordinates(M))
F,V = mergevertices(F,V)

# n = 1
# F,V = subtri(F,V,n; method=:Loop)

s = 0.5
Fs,Vs = scalesimplex(F,V,s)


## Visualisation

strokewidth = 0.1
strokecolor = :white
# Create slightly inward offset version so the visualisation does not co-inside
Vn = V-0.01*vertexnormal(F,V)

fig = Figure(size = (1400,800))

ax1 = Axis3(fig[1, 1], aspect = :data,title="Single scale factor")
mesh!(ax1, GeometryBasics.Mesh(Vn,F), color=:white,transparency=false,shading = FastShading)
hp1 = poly!(ax1, GeometryBasics.Mesh(Vs,Fs), color=:green,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)
# hp1 = mesh!(ax1, GeometryBasics.Mesh(Vs,Fs), color=:green,transparency=false,shading = FastShading)

ax2 = Axis3(fig[1, 2], aspect = :data,title="Per vertex scale factor")
mesh!(ax2, GeometryBasics.Mesh(Vn,F), color=:white,transparency=false,shading = FastShading)
hp2 = poly!(ax2, GeometryBasics.Mesh(Vs,Fs), color=:green,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)
# hp2 = mesh!(ax1, GeometryBasics.Mesh(Vs,Fs), color=:green,transparency=false,shading = FastShading)

ax3 = Axis3(fig[1, 3], aspect = :data,title="Per face scale factor")
mesh!(ax3, GeometryBasics.Mesh(Vn,F), color=:white,transparency=false,shading = FastShading)
hp3 = poly!(ax3, GeometryBasics.Mesh(Vs,Fs), color=:green,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)
# hp3 = mesh!(ax3, GeometryBasics.Mesh(Vs,Fs), color=:green,transparency=false,shading = FastShading)


stepRange = range(1.0,0.0,50)
hSlider = Slider(fig[2, :], range = stepRange, startvalue = 1.0,linewidth=30)

slidercontrol(hSlider,fig)

on(hSlider.value) do s
    X = [v[1] for v in V]
    X = X.-minimum(X)
    X = X./maximum(X)

    VF = simplexcenter(F,V)
    XF = [v[1] for v in VF]
    XF = XF.-minimum(XF)
    XF = XF./maximum(XF)
    
    s1 = s
    s2 = s.*X
    s3 = s.*XF
    Fs1,Vs1 = scalesimplex(F,V,s1)
    Fs2,Vs2 = scalesimplex(F,V,s2)
    Fs3,Vs3 = scalesimplex(F,V,s3)
    
    hp1[1] = GeometryBasics.Mesh(Vs1,Fs1)
    hp2[1] = GeometryBasics.Mesh(Vs2,Fs2)
    hp3[1] = GeometryBasics.Mesh(Vs3,Fs3)
end

set_close_to!(hSlider,0.5)
# fileName = comododir()*"/assets/temp/surface_mesh_smoothing_anim.mp4"
# slider2anim(fig,hSlider,fileName; backforth=true, duration=2)



fig


