using Comodo
using GLMakie
using Statistics
using GeometryBasics
using FileIO

#=
This demo shows the use of the `meshnormal` function to obtain mesh face normal
directions. The demo shows visualisations for a triangular, quadrilateral, and 
a pentagonal mesh. 
=#

fig = Figure(size=(1600,800))

for q=1:1:4
    if q==1
        M=icosahedron()
        F = faces(M)
        V = coordinates(M)
        titleString="triangles"        
    elseif q==2
        M=cube()
        F = faces(M)
        V = coordinates(M)
        titleString="quadrilaterals"        
    elseif q==3
        M=dodecahedron()
        F = faces(M)
        V = coordinates(M)
        titleString="pentagons"        
    elseif q==4
        # Loading a mesh
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)

        # Obtain mesh faces and vertices
        F = togeometrybasics_faces(faces(M))
        V = togeometrybasics_points(coordinates(M))        
        F,V,_ = mergevertices(F,V)    
        M = GeometryBasics.Mesh(V,F)
        titleString="bunny"         
    end

    # Compute mesh vertex normals
    NV = vertexnormal(F,V; weighting="area")
    
    ax1 = Axis3(fig[1, q], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString)
    hp1 = poly!(ax1,M, strokewidth=1,shading=FastShading,color=:white, transparency=false, overdraw=false)    
    # hpa = normalplot(ax1,M,color=:red,linewidth=2)
    # hpa = arrows!(ax1,V,NV,color=:blue)
    hpa = dirplot(ax1,V,NV.*mean(edgelengths(F,V)); color=:blue,linewidth=3,scaleval=1.0,style="from")
end

fig
