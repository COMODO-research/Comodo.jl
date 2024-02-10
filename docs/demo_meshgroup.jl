using Gibbon
using GLMakie
using GeometryBasics
using Statistics
using Rotations
using LinearAlgebra
using FileIO
using SparseArrays

# Example geometry
testCase = 1

r = 1.0
F,V = geosphere(2,r)
for i = 0:1:2
    for j = 0:1:2
        for k = 0:1:2
            if i+j+k>0
                F2,V2 = geosphere(rand(0:2,1)[1],r)
                F2 = map(f-> f.+length(V),F2)
                V2 = map(v-> Point{3, Float64}(i*3*r+v[1],j*3*r+v[2],k*3*r+v[3]),V2)
                append!(F,F2)
                append!(V,V2)
            end
        end
    end
end

C = meshgroup(F)
numGroups = maximum(C)

c = cgrad(:Spectral,numGroups,categorical = true)

M=GeometryBasics.Mesh(V,F)
fig = Figure(size=(800,800))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A multi-object mesh")
poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=2,color=:white, strokecolor=:black, shading = FastShading, transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Grouping")
for q = 1:1:maximum(C)
    poly!(ax2,GeometryBasics.Mesh(V,F[C.==q]), strokewidth=2,color=c[q], strokecolor=:black, shading = FastShading, transparency=false)
end
fig
