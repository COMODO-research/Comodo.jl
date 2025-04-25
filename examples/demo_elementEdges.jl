using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

testCase = 2

if testCase == 1    
    r = 1.0
    n = 1
    Fn1,Vn1 = quadsphere(n,r)
elseif testCase == 2
    r = 1.0
    n = 1
    Fb,Vb = geosphere(n,r)
    E,V,CE,Fb,Cb = tetgenmesh(Fb,Vb)   
end

edgeSet = elementEdges(E)

F = element2faces(E) # Triangular faces

## Visualization
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Mesh edges")
hp1 = mesh!(ax1,GeometryBasics.Mesh(V,F), color=:white, shading = FastShading, transparency=true)
hp2 = wireframe!(ax1,GeometryBasics.Mesh(V,edgeSet), color=:black, linewidth=3)
fig