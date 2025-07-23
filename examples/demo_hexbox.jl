using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics

#=
This demo shows the use of `hexbox` to generate a hexahedral mesh for a 3D box
domain. 
=#

pointSpacing = 0.5
boxDim = [2.5,3.1,4] # Dimensionsions for the box in each direction
boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 

E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)

# Visualisation
cmap = Makie.Categorical(:Spectral) 

Fbs,Vbs = separate_vertices(Fb,V)
Cbs_V = simplex2vertexdata(Fbs,CFb_type)
M = GeometryBasics.Mesh(Vbs,Fbs)

fig = Figure(size=(1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Boundary faces with boundary markers for the hexahedral mesh")
hp2 = meshplot!(ax1, Fbs, Vbs; strokewidth=3, color=Cbs_V, colormap=cmap)
# hp3 = normalplot(ax2,M_Fb; type_flag=:face, color=:black,linewidth=3)

Colorbar(fig[1, 2], hp2)

ax2 = AxisGeom(fig[1, 3], title = "Cut view of hexahedral mesh")
hp3 = meshplot!(ax2, Fbs, Vbs; strokewidth=3, color=:white)

VE  = simplexcenter(E,V)
ZE = [v[3] for v in VE]
Z = [v[3] for v in V]
zMax = maximum(Z)
zMin = minimum(Z)
numSlicerSteps = 3*ceil(Int,(zMax-zMin)/mean(edgelengths(F,V)))

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
slidercontrol(hSlider,ax2)

fig