using Gibbon
using GLMakie
using GeometryBasics
using Statistics
using Rotations

r = 1.0 # radius
n = 2 # Number of refinement iterations
F,V = geoSphere(n,r)

cutLevel = 0.5; 

plateDim = (2,2)
plateElem = (1,1)

FG,VG = quadPlate(plateDim,plateElem)

VC = simplexMiddle(F,V)

L = [mean(V[f])[3]<cutLevel for f ∈ F]

# Define a rotation tensor using Euler angles
Q = RotXYZ(0.0,0.25*pi,0.25*pi)


N,VN=meshnormal(F,V)

## Visualize mesh
fig = Figure(size=(800,800))

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A geodesic sphere n=0")
hp1=poly!(ax1,V,F[.~L], strokewidth=2,color=:white, strokecolor=:black, shading = FastShading, transparency=true)
hp2=poly!(ax1,V,F[L], strokewidth=2,color=:black, strokecolor=:white, shading = FastShading, transparency=true)

hp3=poly!(ax1,GeometryBasics.Mesh(VG,FG), strokewidth=2,color=:red, shading = FastShading, transparency=true)
fig
