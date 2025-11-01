using Comodo
using Comodo.GeometryBasics

f1 = [1,2,3,4,1,5,6,7,8]
f2 = TriangleFace{Int}(1,2,3)
f3 = TriangleFace{Int}(1,2,2)
f4 = QuadFace{Int}(1,2,2,4)
F5 = [TriangleFace{Int}(1,2,3), QuadFace{Int}(1,2,2,4)]

b1 = isunique(f1)
b2 = isunique(f2)
b3 = isunique(f3)
b4 = isunique(f4)
B5 = isunique.(F5)

println("b1=$b1, due to repeated value")
println("b2=$b2, due to the triangle's unique set of node indices")
println("b3=$b3, due to the triangle's non-unique node indices")
println("b4=$b4, due to the quad's non-unique node indices")
println("B5=$B5, since the first entry has unique entries, the second does not")