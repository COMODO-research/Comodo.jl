using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.Statistics
using Comodo.LinearAlgebra
using FileIO

#=
This demo is for the iterative closest point `icp` function which can be used to
align two surfaces. In this demo the surfaces are first "pre-aligned" using 
`surface_align_svd`. This is to help ensure the ICP converges on an aligned 
state which may not happen if the surface is severely misaligned. 
=#

# Get coordinates of set 1
fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
M = load(fileName_mesh)
F1 = tofaces(faces(M))
V1 = topoints(coordinates(M))
F1, V1,_,_ = mergevertices(F1,V1)

s = norm(maxp(V1)-minp(V1)) # rough object size metric
pm1 = surface_centroid(F1, V1)
pc1 = Point{3, Float64}(0.2*s, -0.1*s, 0.25*s)
V1 = [v-pm1 + pc1 for v ∈ V1] 

# Get coordinates of set 2 
F2 = deepcopy(F1)
V2 = deepcopy(V1)

pm2 = surface_centroid(F2, V2)
X = [v[1]-pm2[1] for v in V2]
B = [mean(X[collect(f)])>0.0 for f in F2]
F2, V2 = subtri_dual_local(F2, V2, B)

Q = RotXYZ(0.0,-0.6*π,0.8*π) # Define a rotation tensor using Euler angles

pm2 = surface_centroid(F2, V2)
pc2 = Point{3, Float64}(-0.2*s, 0.1*s, -0.25*s)
V2 = [Point{3, Float64}(Q*(v-pm2)) + pc2 for v ∈ V2] 

# Pre-align using SVD directions (principal component directions)
V2_p1 = surface_align_svd(F1, V1, F2, V2)

# Use ICP to align surface 2 with surface 1
s = 0 # Use smallest (0=default), first (1), or second (2) surface to compote distances to 
n = 100 # Maximum number of iterations
tol = 1e-6 # Tollerance on summed distance difference wrt previous iteration
V2_p2 = icp(V1, V2_p1; n=100, s=0, tol=tol)

# Visualization
GLMakie.closeall()

fig = Figure(size = (1800,1200))
ax1 = AxisGeom(fig[1, 1], title = "Mismatched surfaces")
hp1 = meshplot!(ax1, F1, V1, color=:green, strokecolor=:black, strokewidth=1.0)
hp2 = meshplot!(ax1, F2, V2, color=:red, strokecolor=:black, strokewidth=0.5)

ax2 = AxisGeom(fig[1, 2], title = "SVD alligned surfaces")
hp3 = meshplot!(ax2, F1, V1, color=(:green, 0.75), transparency=true, strokewidth=0.0)
hp4 = meshplot!(ax2, F2, V2_p1, color=(:red, 0.25), transparency=true, strokecolor=:red, strokewidth=1.0)

ax3 = AxisGeom(fig[1, 3], title = "ICP alligned surfaces")
hp5 = meshplot!(ax3, F1, V1, color=(:green, 0.75), transparency=true, strokewidth=0.0)
hp6 = meshplot!(ax3, F2, V2_p2, color=(:red, 0.25), transparency=true, strokecolor=:red, strokewidth=1.0)

Legend(fig[1, 4],[hp1, hp2],["Surface 1", "Surface 2"])

fig