using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using FileIO
using Comodo.LinearAlgebra

GLMakie.closeall()

## Loading example mesh
fileName_mesh = joinpath(comododir(),"assets","stl","david.stl")
M = load(fileName_mesh)

# Obtain mesh faces and vertices
F = tofaces(faces(M))
V = topoints(coordinates(M))
F,V,_ = mergevertices(F,V)

## Seeding points 

# Define start point, here top point is used
z = [v[3] for v ∈ V]
indStart = [findmax(z)[2]]

# Number of desired points
numPoints = 250  

# Do distance marching based point seeding 
ind,d,l = distseedpoints(F, V, numPoints)   

# Convert seeds to downsampled mesh 
Fp, Vp = seedpoints2mesh(F, V, ind, l)

## visualisation
cMap_dist = :bluesreds
cMap_cat = :Spectral #Makie.Categorical(:Spectral) 
fig = Figure(size=(1800, 600))

ax1 = AxisGeom(fig[1, 1], title = "Distances")
hp1 = meshplot!(ax1, F, V, color=d, colormap=cMap_dist, strokewidth=0.0)
hp2 = scatter!(ax1, V[ind], color=:black, markersize=10, depth_shift=-0.01f0)
Colorbar(fig[1, 2], hp1)

ax2 = AxisGeom(fig[1, 3], title = "Point regions, n=$numPoints points")
hp3 = meshplot!(ax2, F, V, color=l, colormap=cMap_cat, strokewidth = 0.1)
hp4 = scatter!(ax2, V[ind], color=:black, markersize=10, depth_shift=-0.01f0)
Colorbar(fig[1, 4], hp3)

ax3 = AxisGeom(fig[1, 5], title = "Reconstructed surface")
Fps, Vps = separate_vertices(Fp, Vp)
hp5 = meshplot!(ax3, Fps, Vps, strokewidth=0.5)
# hp6 = scatter!(ax3, V[ind], color=:black, markersize=15, depth_shift=-0.01f0)

stepRange = 1:1:1000
hSlider = Slider(fig[2, :], range = stepRange, startvalue = numPoints, linewidth=30)

on(hSlider.value) do numPoints 
    ind,d,l = distseedpoints(F,V,numPoints; ind=[1])   

    Fp, Vp = seedpoints2mesh(F, V, ind, d, l)
    
    ax2.title = "Point regions, n=$numPoints points"
    hp1.color = d
    hp2[1] = V[ind]

    hp3.color = l
    hp4[1] = V[ind]

    Fps, Vps = separate_vertices(Fp, Vp)
    hp5[1] = GeometryBasics.Mesh(Vps, Fps)
    # hp6[1] = V[ind]
end

slidercontrol(hSlider,ax1)

screen = display(GLMakie.Screen(), fig)


