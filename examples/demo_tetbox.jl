using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics

#=
This demo shows the use of `tetbox` to generate a tetrahedral mesh for a 3D box
domain. 
=#

boxDim = [2.5,3.1,4] # Dimensions for the box in each direction
pointSpacing = 0.5

E, V, Fb, Cb = tetbox(boxDim,pointSpacing)

# Visualisation
GLMakie.closeall()

cmap = Makie.Categorical(:Spectral) 

Fbs,Vbs = separate_vertices(Fb,V)
Cbs_V = simplex2vertexdata(Fbs,Cb)
M = GeometryBasics.Mesh(Vbs,Fbs)

fig = Figure(size=(1600,800))

ax1 = AxisGeom(fig[1, 1][1, 1], title = "Boundary faces with boundary markers for the tetrahedral mesh")
hp2 = meshplot!(ax1, Fbs, Vbs, strokewidth=2, color=Cbs_V, colorrange = (1,6), colormap=cmap)
# hp3 = normalplot(ax1,Fb,V; type_flag=:face, color=:black,linewidth=3)
Colorbar(fig[1, 1][1, 2], hp2)

ax2 = AxisGeom(fig[1, 2], title = "Cut view of tetrahedral mesh")
hp3 = meshplot!(ax2, Fbs, Vbs, strokewidth=2, color=:white)

VE  = simplexcenter(E,V)
ZE = [v[3] for v in VE]
Z = [v[3] for v in V]
zMax = maximum(Z)
zMin = minimum(Z)
numSlicerSteps = 3*ceil(Int,(zMax-zMin)/mean(edgelengths(Fb,V)))

stepRange = range(zMin,zMax,numSlicerSteps)
hSlider = Slider(fig[2, :], range = stepRange, startvalue = mean(stepRange),linewidth=30)

on(hSlider.value) do z 

    B = ZE .<= z
    indShow = findall(B)
    if isempty(indShow)
        hp3.visible=false        
    else        
        hp3.visible=true
        Fs = element2faces(E[indShow])
        Fs,Vs = separate_vertices(Fs,V)
        Ms = GeometryBasics.Mesh(Vs,Fs)
        hp3[1] = Ms
    end

end
# hSlider.selected_index[]+=1
slidercontrol(hSlider,ax2)

fig
