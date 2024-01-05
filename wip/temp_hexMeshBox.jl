using Gibbon
using GLMakie
using GeometryBasics


boxDim = [2,2,2]
boxEl = [3,3,3]

E,V,F,Fb,CFb_type = hexMeshBox(boxDim,boxEl)   

M = GeometryBasics.Mesh(V,F)
Mb = GeometryBasics.Mesh(V,Fb)

# N,VN=meshnormal(M)

#Visualize mesh
fig = Figure(size = (800,800))

ax=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Hexahedral mesh")
hp=poly!(ax,Mb,strokewidth=5,color=:red, transparency=true, overdraw=false)
# hp=poly!(ax,M,strokewidth=5,color=CF_type, transparency=false, overdraw=false,colormap = Reverse(:Spectral), colorrange=(0,6))

# hp=poly!(ax,[GeometryBasics.Mesh(VM,f) for f ∈ eachrow(Fb)],strokewidth=5,color=CFb_type, transparency=true, overdraw=false,colormap = Reverse(:Spectral), colorrange=(0,6))

# for i ∈ eachindex(F)    
#    poly!(ax,GeometryBasics.Mesh(VM,F[i]),strokewidth=5,color=CF_type[i], transparency=false, overdraw=false,colormap = Reverse(:Spectral), colorrange=(0,6))
# end

# for i ∈ eachindex(Fb)    
#    poly!(ax,GeometryBasics.Mesh(V,Fb[i]),strokewidth=5,color=CFb_type[i], transparency=false, overdraw=false,colormap = Reverse(:Spectral), colorrange=(0,6))
# end

# Colorbar(fig[1, 1]) 
# hpa=arrows!(ax,VN,N)
fig
