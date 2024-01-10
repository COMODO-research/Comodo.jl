using Gibbon
using GLMakie
using Meshes
# using GeometryBasics

points = Meshes.Point3[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
connec = Meshes.connect.([(1,4,3,2),(5,6,7,8),(1,2,6,5),(3,4,8,7),(1,5,8,4),(2,3,7,6)])
mesh   = SimpleMesh(points, connec)

function refdn(mesh,n)
    for _ = 1:1:n
        mesh = refine(mesh, CatmullClark())
    end
    return mesh
end

n=5
ref = refdn(mesh,n)

fig = Figure(size = (1600, 800))
viz(fig[1,1], ref, showfacets = true)

fig

