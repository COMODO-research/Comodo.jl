using Gibbon
using TetGen
using GeometryBasics
using GeometryBasics: Mesh, QuadFace
using GLMakie

# function tetgenResult2Mesh(R)

# end

# Construct a cube out of Quads
points = Point{3, Float64}[
    (0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
    (2.0, 2.0, 0.0), (0.0, 2.0, 0.0),
    (0.0, 0.0, 2.0), (2.0, 0.0, 2.0),
    (2.0, 2.0, 2.0), (0.0, 2.0, 2.0)
]

facets = QuadFace{Cint}[
    1:4,
    5:8,
    [1,5,6,2],
    [2,6,7,3],
    [3, 7, 8, 4],
    [4, 8, 5, 1]
]

markers = Cint[-1, -2, 0, 0, 0, 0]
# attach some additional information to our faces!
mesh = Mesh(points, meta(facets, markers=markers))
result = tetrahedralize(mesh, "vpq1.414a0.1")

function tetgen2facemesh(result)
    E = faces(result) # elements
    V = coordinates(result) # Vertices

    # Get faces
    F = Vector{TriangleFace{Int64}}(undef,length(E)*4) # Initialise faces
    for q = eachindex(E) # Loop over all elements
        e = E[q] # The current element 
        qf = 1 + (q-1)*4 # Index mapping into face array 
        F[qf  ] = TriangleFace{Int64}(e[1],e[2],e[3]) 
        F[qf+1] = TriangleFace{Int64}(e[1],e[2],e[4]) 
        F[qf+2] = TriangleFace{Int64}(e[2],e[3],e[4])
        F[qf+3] = TriangleFace{Int64}(e[3],e[1],e[4])
    end

    return GeometryBasics.Mesh(V,F) # Return GeometryBasics mesh type
end

M = tetgen2facemesh(result)

# N,VN=meshnormal(M)

#Visualize mesh
fig = Figure(size=(1200,800))

ax2=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Domain to mesh")
poly!(ax2,GeometryBasics.Mesh(points,facets),strokewidth=5, color=:blue, transparency=true, overdraw=false,shading = NoShading)
# arrows!(ax2,VN,N)

ax2=Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Tetgen tetrahedral mesh")
poly!(ax2,M,strokewidth=5, color=:gray, transparency=false, overdraw=false,shading = NoShading)
# arrows!(ax2,VN,N)


fig

# #Visualize mesh
# GLMakie.activate!(inline=false) # To avoid plotting in plotpane as per: https://github.com/MakieOrg/Makie.jl/issues/2956
# fig = Figure()
# ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Tetrahedron mesh")

# mesh!(ax1,normal_mesh(result), color=:blue, transparency=true, overdraw=false)
# wireframe!(ax1,result)

# fig