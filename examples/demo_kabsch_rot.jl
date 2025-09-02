using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using FileIO

#=
This demo is for the `kabsch_rot` function. The Kabsch method allows one to 
determine the rotation that occurred between two meshes (or point sets) with 
point-to-point correspondence. In this demo a triangulated surface is loaded 
from an STL file. Next the coordinates are rotated, and the Kabsch method is 
used to retrieve the rotation performed, allowing one to "unrotate" the data.  
=#

# Loading a mesh

fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
M = load(fileName_mesh)

# Obtain mesh faces and vertices
F = tofaces(faces(M))
V1 = topoints(coordinates(M))
F,V1,_,_ = mergevertices(F,V1)

# Define a rotation tensor using Euler angles
Q = RotXYZ(0.0,0.25*π,0.25*π)
V2 = [GeometryBasics.Point{3, Float64}(Q*v) for v ∈ V1] 

R = kabsch_rot(V2,V1)
V3 = [GeometryBasics.Point{3, Float64}(R*v) for v ∈ V2] 

# Visualization
GLMakie.closeall()

fig = Figure(size = (800,800))
ax = AxisGeom(fig[1, 1])
hp1 = meshplot!(ax, F, V1, color=:green)
hp2 = meshplot!(ax, F, V2, color=:red)
hp3 = edgeplot!(ax, F, V3, color=:red, linewidth=1)
Legend(fig[1, 2],[hp1,hp2,hp3],["Initial","Rotated","Back rotated"])
fig