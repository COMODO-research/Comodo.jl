using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.GLMakie.FileIO

#=
This demo shows the use of the `meshnormal` function to obtain mesh face normal
directions. The demo shows visualisations for a triangular, quadrilateral, and 
a pentagonal mesh. 
=#

GLMakie.closeall()

fig = Figure(size=(1600,800))

for q=1:1:4
    if q==1
        F,V = icosahedron()
        titleString="triangles"        
    elseif q==2
        F,V = cube()
        titleString="quadrilaterals"        
    elseif q==3
        F,V = dodecahedron()        
        titleString="pentagons"        
    elseif q==4
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = tofaces(faces(M))
        V = topoints(coordinates(M))        
        F,V,_ = mergevertices(F,V)    
        M = GeometryBasics.Mesh(V,F)
        titleString="bunny"         
    end

    # Compute mesh vertex normals
    NV = vertexnormal(F,V; weighting=:size)
    
    ax1 = AxisGeom(fig[1, q], title = titleString)
    hp1 = meshplot!(ax1, F, V)    
    # hpa = normalplot(ax1,M,color=:red,linewidth=2)    
    hpa = dirplot(ax1,V,NV.*mean(edgelengths(F,V)); color=:blue,linewidth=3,scaleval=1.0,style=:from)
end

fig