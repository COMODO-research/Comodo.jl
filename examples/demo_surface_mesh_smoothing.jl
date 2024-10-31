using Comodo, GLMakie, GeometryBasics, FileIO, Statistics, Random

Random.seed!(1) # Set seed so demo performs the same each time

testCase = 1 
if testCase == 1 # Triangle mesh bunny
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)
    F = faces(M)
    V = coordinates(M)    
    V = topoints(V) 
    F = tofaces(F)
    F,V,_,_ = mergevertices(F,V)
elseif testCase == 2 # Quad mesh sphere
    F,V = quadsphere(4,100)
end

V0 = deepcopy(V)

V = map(x-> x.+ (5.0 .* randn(eltype(V))),V)

# Set smoothing parameters
tol = 1e-3
nMax = 15 # Maximum number of iterations

λ = 0.5 # Laplacian smoothing parameter

# HC parameters
α = 0.1
β = 0.5

# Use constrained nodes for the bottom half
testMode = 2
if testMode == 2
    z = [v[3] for v in V]
    constrained_points = findall(z.<mean(z))
else
    constrained_points = nothing
end

## VISUALISATION

strokeWidth1 = 1

Vs_HC = smoothmesh_hc(F,V,nMax,α,β; tolDist=tol, constrained_points = constrained_points)
Vs_LAP = smoothmesh_laplacian(F,V,nMax,λ; constrained_points = constrained_points)
Ds = [sqrt(sum((Vs_HC[i]-V0[i]).^2)) for i ∈ eachindex(V)]
Ds_LAP = [sqrt(sum((Vs_LAP[i]-V0[i]).^2)) for i ∈ eachindex(V)]
cLim = maximum(Ds_LAP).*(0.0,1.0)

fig = Figure(size=(1200,1200))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original")
hp1=poly!(ax1,GeometryBasics.Mesh(V0,F),strokewidth=strokeWidth1,color=:white, shading = FastShading)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Noisy")
hp2=poly!(ax2,GeometryBasics.Mesh(V,F),strokewidth=strokeWidth1,color=:white, shading = FastShading)

ax3 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")
hp3 = poly!(ax3,GeometryBasics.Mesh(Vs_LAP,F),strokewidth=strokeWidth1,color=Ds_LAP, shading = FastShading,colormap=Makie.Reverse(:Spectral),colorrange=cLim)

ax4 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")
hp4 = poly!(ax4,GeometryBasics.Mesh(Vs_HC,F),strokewidth=strokeWidth1,color=Ds, shading = FastShading,colormap=Makie.Reverse(:Spectral),colorrange=cLim)

Colorbar(fig[:, 3],hp3,label = "Distance")

stepRange = 0:1:nMax
hSlider = Slider(fig[3, :], range = stepRange, startvalue = 0,linewidth=30)

slidercontrol(hSlider,fig)

on(hSlider.value) do stepIndex
    # Update first plot
    Vs_LAP = smoothmesh_laplacian(F,V,stepIndex,λ; constrained_points = constrained_points)
    Ds_LAP = [sqrt(sum((Vs_LAP[i]-V0[i]).^2)) for i ∈ eachindex(V)]
    ax3.title= "Laplacian smoothed n = " * string(stepIndex) * " times"
    hp3.color = Ds_LAP
    hp3[1] = GeometryBasics.Mesh(Vs_LAP,F)

    # Update second plot
    Vs_HC = smoothmesh_hc(F,V,stepIndex,α,β; tolDist=tol, constrained_points = constrained_points)
    Ds_HC = [sqrt(sum((Vs_HC[i]-V0[i]).^2)) for i ∈ eachindex(V)]
    ax4.title= "HC smoothed n = " * string(stepIndex) * " times"
    hp4.color = Ds_HC
    hp4[1] = GeometryBasics.Mesh(Vs_HC,F)
end

set_close_to!(hSlider,0)

# fileName = comododir()*"/assets/temp/surface_mesh_smoothing_anim.mp4"
# slider2anim(fig,hSlider,fileName; backforth=true, duration=2)

fig