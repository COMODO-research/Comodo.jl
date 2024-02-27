using Comodo
using GLMakie
using GeometryBasics
using Statistics
using Rotations
using LinearAlgebra
using FileIO
using SparseArrays

# Example geometry
F,V = geosphere(3,1.0)
VC = simplexcenter(F,V)
F = [F[i] for i in findall(map(v-> v[3]>0,VC))] # Remove some faces

Fc,Vc = remove_unused_vertices(F,V)


M = GeometryBasics.Mesh(Vc,Fc)
fig = Figure(size=(1200,1200))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A sliced mesh")


hp3 = poly!(ax1,M, strokewidth=1,color=:white, strokecolor=:blue, shading = FastShading, transparency=false)
hp3 = normalplot(ax1,M)

fig
