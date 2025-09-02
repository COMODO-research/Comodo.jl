using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.GLMakie.Colors
using FileIO

#=
This demo shows the use of the `facenormal` function to obtain mesh face normal
directions. The demo shows visualisations for a triangular, quadrilateral, and 
a pentagonal mesh. 
=#

GLMakie.closeall()

fig = Figure(size=(1600,800))
cAlpha = RGBA(1.0, 1.0, 1.0,0.25)
for q=1:1:5
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
        F,V = dodecahedron()
        invert_faces!(F)
        titleString="inverted pentagons"   
    elseif q==5             
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh)
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        F,V,_,_ = mergevertices(F,V)
        titleString="Stanford bunny"   
    end

    # Compute mesh face normals
    N = facenormal(F,V)    
    VN = simplexcenter(F,V)
    
    ax1 = AxisGeom(fig[1, q], title = titleString)
    hp1 = meshplot!(ax1,F,V, color=cAlpha, transparency=true)
    # hpa = arrows!(ax1,VN,N,color=:blue)  
    hpa = dirplot(ax1,VN,N; color=:blue,linewidth=3, scaleval=1.0, style=:from)
end

fig