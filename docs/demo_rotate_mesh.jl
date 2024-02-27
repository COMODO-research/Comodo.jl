using Comodo
using GLMakie
using GeometryBasics
using FileIO
using Rotations

"""
In this demo a mesh is loaded, in this case a triangulated surface from an STL 
file. Next the coordinates are rotated, and the unrotated and rotated meshes 
are visualized. 
"""

# Loading a mesh
fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
M = load(fileName_mesh)

# Obtain mesh faces and vertices
F = faces(M)
V = coordinates(M)

# Define a rotation tensor using Euler angles
Q = RotXYZ(0.0,0.25*π,0.25*π)

# Rotate the coordinates
V2 = [GeometryBasics.Point{3, Float64}(Q*v) for v ∈ V] 

# Create a mesh description for visualisation 
Mr = GeometryBasics.Mesh(V2,F); 

fig = Figure(size = (800,800))
ax = Axis3(fig[1, 1], aspect = :data)

hp1 = poly!(ax, M, strokewidth=2,color=:green,transparency=true,shading = FastShading)
hp2 = poly!(ax, Mr, strokewidth=2,color=:red,shading = FastShading)

Legend(fig[1, 2],[hp1,hp2],["Initial","Rotated"])

fig