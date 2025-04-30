using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.LinearAlgebra

#=
This demo shows the use of `hexbox` to generate a hexahedral mesh for a 3D box
domain. 
=#

r = 1.0
f = 0.4
pointSpacing = r/10.0

function hexsphere(r,pointSpacing; f=0.4)    
    boxDim = fill(2*f*r,3) # Dimensionsions for the box in each direction
    boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 
    E,V,_,Fb,_ = hexbox(boxDim,boxEl)
    ind = unique(reduce(vcat,Fb))    
    V2 = [v.* (r/norm(v)) for v in V[ind]]
    numSteps = ceil(Int,((1.0-f)*r)/pointSpacing)    
    append!(E,fromtomesh!(Fb,V,V2,numSteps; correspondence=:faces))
    return E,V
end

E,V = hexsphere(r,pointSpacing; f=0.4)    
Fn = element2faces(E)

# Visualisation
cmap = Makie.Categorical(:Spectral) 
a = 0.25
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Volumetric mesh (cut view)")
hp3 = poly!(ax1, GeometryBasics.Mesh(V,Fn), strokewidth=3,shading=FastShading,strokecolor=:black, color=:white, transparency=false, overdraw=false)

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

