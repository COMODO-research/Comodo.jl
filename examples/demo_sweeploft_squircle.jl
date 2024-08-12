using Comodo
using GeometryBasics
using GLMakie
using Rotations
using Statistics
using LinearAlgebra

function makeMesh(τ = 1.0)

    nc = 150 # Number of points on guide curve
    r = 5
    Vc = [GeometryBasics.Point{3, Float64}(r*cos(t),r*sin(t),0.0) for t ∈ range(0,2*π,nc)]

    # Define section curves
    np = 40 # Number of section points

    # Section 1
    V1 = squircle(2.0,np,τ; dir=:acw)
    # V1 = evenly_sample(V1, np; close_loop=true)
    V1 = [Point{3,Float64}(v[1],v[3],v[2]) for v ∈ V1]
    V1 = [v .+ Vc[end] for v ∈ V1]

    # Section 2
    V2 = squircle(2.0,np,τ; dir=:acw)
    # V2 = evenly_sample(V2, np; close_loop=true)
    V2 = [Point{3,Float64}(v[1],v[3],v[2]) for v ∈ V2]
    V2 = [v .+ Vc[end] for v ∈ V2]

    F,V = sweeploft(Vc,V1,V2; face_type=:quad, num_twist=3)
    F,V = mergevertices(F,V)
    # F = invert_faces(F)   
    # C = ceil.(collect(1:length(V))./(length(V1)))
    K1,K2,U1,U2,H,G = mesh_curvature_polynomial(F,V; growsteps=2)

    return F,V,G
end

F,V,C = makeMesh(1.0)

#########


# Visualization

cMap = :Spectral

markersize = 6
linewidth = 2

fig = Figure(size = (1600,1600))
ax = Axis3(fig[1, 1],aspect = :data,title="Swept lofting")

stepRange1 = range(0.0,1.0,100)
hSlider1 = Slider(fig[2, 1], range = stepRange1, startvalue = 0,linewidth=30)

# hp1 = poly!(ax, GeometryBasics.Mesh(V,F), strokecolor=:black, strokewidth=0.5,color=C,transparency=false,shading = FastShading,colormap=cMap)
hp1 = mesh!(ax, GeometryBasics.Mesh(V,F), color=C,transparency=false,shading = FastShading,colormap=cMap)

on(hSlider1.value) do τ
    F,V,C = makeMesh(τ)
    hp1[1] = GeometryBasics.Mesh(V,F)
    hp1.color = C    
end

fig

# fileName = comododir()*"/assets/img/sweep_loft_squircle_anim_02.mp4"
# slider2anim(fig,hSlider1,fileName; backforth=true, duration=4)