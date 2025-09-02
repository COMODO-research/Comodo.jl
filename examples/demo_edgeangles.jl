using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

F,V = cube(sqrt(3))

# Build deformation gradient tensor to induce shear with known angles
f = zeros(3,3)
for i = 1:3
    f[i,i] = 1.0
end
a = pi/4 # "45 degree shear"  
f[1,2] = tan(a) 

V2 = [eltype(V)(f*v) for v ∈ V] # Shear the cube

A = edgeangles(F,V; deg = true)
A2 = edgeangles(F,V2; deg = true)

## Visualize mesh
GLMakie.closeall()

Fs,Vs = separate_vertices(F,V2) # Separate to enable per face color shading
fig = Figure(size=(800,800))
ax1 = AxisGeom(fig[1, 1], title = "Edge angles")
hp3 = meshplot!(ax1, Fs, Vs, strokewidth=3, color=reduce(vcat,A2), colormap=Makie.Reverse(:Spectral))
Colorbar(fig[1,2],hp3)
fig
