using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 2.0 # radius
n = 3 # Number of refinement iterations
F,V = geosphere(n,r)
pointSpacing = pointspacingmean(F,V)

n2 = n+2
pointSpacingDesired = pointSpacing/4.0

println("pointSpacingDesired ",pointSpacingDesired)
NV = spacing2numvertices(F,V,pointSpacingDesired)
println("num V ",NV)

F2,V2 = geosphere(n2,r)
println("num true ",length(V2))

pointSpacingTrue = pointspacingmean(F2,V2)

println("pointSpacingTrue ",pointSpacingTrue)