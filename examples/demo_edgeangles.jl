using Comodo
using GLMakie
using GeometryBasics

M = cube(1.0)
F = faces(M)
V = coordinates(M)

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
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Edge angles")
hp2 = poly!(ax1,GeometryBasics.Mesh(V2,F), strokewidth=3,color=:white, shading = FastShading)

fig
