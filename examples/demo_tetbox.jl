using Comodo
using GLMakie
using Statistics

#=
This demo shows the use of `tetbox` to generate a tetrahedral mesh for a 3D box
domain. 
=#


boxDim = [2.5,3.1,4] # Dimensions for the box in each direction
pointSpacing = 0.5

E, V, Fb, Cb = tetbox(boxDim,pointSpacing)

# Visualisation
cmap = Makie.Categorical(:Spectral) 

Fbs,Vbs = separate_vertices(Fb,V)
Cbs_V = simplex2vertexdata(Fbs,Cb)
M = GeometryBasics.Mesh(Vbs,Fbs)

fig = Figure(size=(1600,800))

ax1 = Axis3(fig[1, 1][1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary faces with boundary markers for the tetrahedral mesh")
hp2 = poly!(ax1,M, strokewidth=3,shading=FastShading,strokecolor=:black, color=Cbs_V, transparency=false, overdraw=false,colorrange = (1,6),colormap=cmap)
# hp3 = normalplot(ax1,Fb,V; type_flag=:face, color=:black,linewidth=3)

ax2 = Axis3(fig[1, 1][1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Cut view of tetrahedral mesh")
hp3 = poly!(ax2,M, strokewidth=3,shading=FastShading,strokecolor=:black, color=:white, transparency=false, overdraw=false)

VE  = simplexcenter(E,V)
ZE = [v[3] for v in VE]
Z = [v[3] for v in V]
zMax = maximum(Z)
zMin = minimum(Z)
numSlicerSteps = 3*ceil(Int,(zMax-zMin)/mean(edgelengths(Fb,V)))

stepRange = range(zMin,zMax,numSlicerSteps)
hSlider = Slider(fig[2, 1], range = stepRange, startvalue = mean(stepRange),linewidth=30)

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
