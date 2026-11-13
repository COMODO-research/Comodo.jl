using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.LinearAlgebra

#=
This demo shows the use of `hexhemisphere` to generate a hexahedral mesh for a 
hemisphere.  
=#

##

nSmooth = 50
r = 5.0
f = 0.75
pointSpacing = 0.75
nSmooth = 5

E, V, F, Fb, Cb = hexhemisphere(r,pointSpacing; f=f, nSmooth=nSmooth)    

##
Fn = element2faces(E)

# Visualisation
GLMakie.closeall()

cmap = Makie.Categorical(Reverse(:Spectral)) 

fig = Figure(size=(800,800))

Fbs,Vbs = separate_vertices(Fb,V)
Cbs_V = simplex2vertexdata(Fbs, Cb, Vbs)

ax1 = AxisGeom(fig[1, 1],  title = "Boundary surface")
hp1 = meshplot!(ax1, Fbs, Vbs, strokewidth=1, color=Cbs_V, colormap=cmap)
Colorbar(fig[1, 2], hp1)

ax2 = AxisGeom(fig[1, 3],  title = "Volumetric mesh (cut view)")
hp3 = meshplot!(ax2, Fbs, Vbs, strokewidth=1, color=:white)


VE  = simplexcenter(E,V)
ZE = [v[3] for v in VE]
Z = [v[3] for v in V]
zMax = maximum(Z)
zMin = minimum(Z)
numSlicerSteps = 3*ceil(Int,(zMax-zMin)/pointSpacing)

stepRange = range(zMin,zMax,numSlicerSteps)
hSlider = Slider(fig[2, :], range = stepRange, startvalue = stepRange[end],linewidth=30)

on(hSlider.value) do z     
    indShow = findall(ZE .<= z)
    if isempty(indShow)
        hp3.visible=false        
    else        
        hp3.visible=true
        Fn = element2faces(E[indShow])        
        Fns,Vns = separate_vertices(Fn,V)
        hp3[1] = GeometryBasics.Mesh(Vns,Fns)
    end
end

slidercontrol(hSlider,ax1)

fig

