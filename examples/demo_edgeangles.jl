using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

F,V = cube(1.0)

# Build deformation gradient tensor to induce shear with known angles
f = zeros(3,3)
for i=1:3
    f[i,i]=1.0
end
a = pi/4 # "45 degree shear"  
f[1,2] = tan(a) 

V2 = [eltype(V)(f*v) for v ∈ V] # Shear the cube

A = edgeangles(F,V)
A2 = edgeangles(F,V2)


## Visualize mesh

Fs,Vs = separate_vertices(F,V2)

fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Edge angles")
hp3 = poly!(ax1,GeometryBasics.Mesh(Vs,Fs), strokewidth=3,color=reduce(vcat,A2), shading = FastShading,colormap=Makie.Reverse(:Spectral))
Colorbar(fig[1,2],hp3)
fig
