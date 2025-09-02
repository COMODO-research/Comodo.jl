using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.LinearAlgebra

#=
This demo shows the use of `hexsphere` to generate a hexahedral mesh for a 3D 
sphere domain. 
=#

nSmooth = 50
r = 5.0
f = 0.75
pointSpacing = 0.5
E,V = hexsphere(r,pointSpacing; f=f, nSmooth=nSmooth)    
Fn = element2faces(E)

# Visualisation
GLMakie.closeall()

cmap = Makie.Categorical(:Spectral) 
a = 0.25
fig = Figure(size=(800,800))

ax1 = AxisGeom(fig[1, 1],  title = "Volumetric mesh (cut view)")
hp3 = meshplot!(ax1, Fn, V, strokewidth=2)

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

