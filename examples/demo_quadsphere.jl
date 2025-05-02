using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Statistics

#=
This demo shows the use of `quadsphere` to generate a quadrilateral surface mesh
for a 3D sphere domain. 
=#

r = 1.0 # Radius

pointSpacing = 0.5
F,V,C = quadsphere(r,pointSpacing)

## Visualization
strokewidth = 2
cmap = cgrad(:Spectral,6, categorical = true)
Fs,Vs = separate_vertices(F,V)
Cs = simplex2vertexdata(Fs,C)

fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "quadrangulated sphere, point spacing = $pointSpacing")
hp1 = poly!(ax1,GeometryBasics.Mesh(Vs,Fs), strokewidth=strokewidth,color=Cs,
shading=FastShading,transparency=false, colorrange = (1,6),colormap=cmap)
Colorbar(fig[1, 2], hp1)

stepRange = (2.0*π*r)/4.0 ./ collect(1:1:30)
hSlider = Slider(fig[2, 1], range = stepRange, startvalue = pointSpacing,linewidth=30)

on(hSlider.value) do pointSpacing
    F,V,C = quadsphere(r,pointSpacing)
    Fs,Vs = separate_vertices(F,V)
    Cs = simplex2vertexdata(Fs,C)

    hp1[1] = GeometryBasics.Mesh(Vs,Fs)
    hp1.color = Cs
    ax1.title = "quadrangulated sphere, point spacing = $pointSpacing"
end
slidercontrol(hSlider,ax1)

fig