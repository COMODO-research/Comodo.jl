using Gibbon
using GeometryBasics
using BenchmarkTools
using SparseArrays
using GLMakie

r = 1 
M = platonicsolid(4,r) # Get icosahedron
V = coordinates(M) # Get vertices
F = faces(M) # Get faces
F = F[1:end-1]

# F,V = subTri(F,V,2)

E = meshEdges(F)

# FACE-EDGE connectivity
# con_f2e = [indReverse[[1,2,3].+ (i-1)*3] for i ∈ eachindex(F)]


E_uni,indUni,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)

con_F2E = [Vector{Int64}(a) for a ∈ eachrow(reshape(indReverse,3,length(F)))] # [indReverse[[1,2,3].+ (i-1)*3] for i ∈ eachindex(F)]
# con_F2E = reshape(indReverse,3,length(F)) # [indReverse[[1,2,3].+ (i-1)*3] for i ∈ eachindex(F)]
# ind = [i.*ones(Int64,3) for i=1:length(F)]
ind = repeat(1:1:length(F),3); 
A = sparse(indReverse,ind,ind)

con_E2F = Vector{Vector{Int64}}(undef,length(E_uni))
for q ∈ eachindex(E_uni)
    _,v = findnz(A[q,:])
    con_e2f[q] = v
end

# con_e2f[indReverse]

fig = Figure(size=(1600,800))

ax=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")
hp=poly!(ax,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white, shading=FastShading, overdraw=false)
indE=21
hp=lines!(ax,[V[E_uni[indE][1]],V[E_uni[indE][2]]], color=:red,linewidth = 10)
hp=poly!(ax,GeometryBasics.Mesh(V,F[9]), strokewidth=4,color=:blue, shading=FastShading, overdraw=false)

fig

