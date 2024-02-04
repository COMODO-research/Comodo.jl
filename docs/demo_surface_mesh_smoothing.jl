using Gibbon, GLMakie, GeometryBasics, FileIO, Statistics, Random

Random.seed!(1) # Set seed so demo performs the same each time

# Loading a mesh
fileName_mesh = joinpath(gibbondir(),"assets","stl","stanford_bunny_low.stl")
M = load(fileName_mesh)
F = faces(M)
V = coordinates(M)
F = togeometrybasics_faces(F) 
V = togeometrybasics_points(V) 
F,V,ind1,ind2 = mergevertices(F,V; roundVertices=true)

V0 = deepcopy(V)

M0 = GeometryBasics.Mesh(V0,F)
V = map(x-> x.+ (5.0 .* randn(eltype(V))),V)
M = GeometryBasics.Mesh(V,F)

E = meshedges(F)
E_uni,_,_ = unique_simplices(E)
con_V2V = con_vertex_vertex(E_uni)

# Set smoothing parameters

# Laplacian smoothing parameters
λ = 0.5

# HC parameters
α = 0.1
β = 0.5
tol = 1e-3

nMax = 100 # Maximum number of iterations


## VISUALISATION

strokeWidth1 = 1

Vs_HC = smoothmesh_hc(F,V,con_V2V; n=nMax,α=α,β=β,tolDist=tol)
Vs_LAP = smoothmesh_laplacian(F,V,con_V2V; n=nMax,λ=λ)
Ds = [sqrt(sum((Vs_HC[i]-V0[i]).^2)) for i ∈ eachindex(V)]
Ds_LAP = [sqrt(sum((Vs_LAP[i]-V0[i]).^2)) for i ∈ eachindex(V)]
cLim = maximum(Ds_LAP).*(0.0,1.0)

fig = Figure(size=(900,900))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original")
hp1=poly!(ax1,M0,strokewidth=strokeWidth1,color=:white, shading = FastShading)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Noisy")
hp2=poly!(ax2,M,strokewidth=strokeWidth1,color=:white, shading = FastShading)

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
    Vs_LAP = smoothmesh_laplacian(F,V,con_V2V; n=stepIndex,λ=λ)
    Ds_LAP = [sqrt(sum((Vs_LAP[i]-V0[i]).^2)) for i ∈ eachindex(V)]
    ax3.title= "Laplacian smoothed n = " * string(stepIndex) * " times"
    hp3.color = Ds_LAP
    hp3[1] = GeometryBasics.Mesh(Vs_LAP,F)

    # Update second plot
    Vs_HC = smoothmesh_hc(F,V,con_V2V; n=stepIndex,α=α,β=β,tolDist=tol)
    Ds_HC = [sqrt(sum((Vs_HC[i]-V0[i]).^2)) for i ∈ eachindex(V)]
    ax4.title= "HC smoothed n = " * string(stepIndex) * " times"
    hp4.color = Ds_HC
    hp4[1] = GeometryBasics.Mesh(Vs_HC,F)
end

set_close_to!(hSlider,0)

fig