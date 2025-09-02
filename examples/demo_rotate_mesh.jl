using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using FileIO

#=
In this demo a mesh is loaded, in this case a triangulated surface from an STL 
file. Next the coordinates are rotated, and the unrotated and rotated meshes 
are visualized. 
=#

# Loading a mesh
fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
M = load(fileName_mesh)

# Obtain mesh faces and vertices
F = [TriangleFace{Int}(f) for f in faces(M)]
V = coordinates(M)

# Define a rotation tensor using Euler angles
Q = RotXYZ(0.0,0.25*π,0.25*π)

# Rotate the coordinates
V2 = [GeometryBasics.Point{3, Float64}(Q*v) for v ∈ V] 

# Visualisation
GLMakie.closeall()

fig = Figure(size = (800,800))
ax = AxisGeom(fig[1, 1])
hp1 = meshplot!(ax, F,  V, color=:green)
hp2 = meshplot!(ax, F, V2, color=:red)
Legend(fig[1, 2],[hp1,hp2],["Initial","Rotated"])
fig