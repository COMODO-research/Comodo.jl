using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.GLMakie.Colors

#=
This demo shows the use of the `facenormal` function to obtain mesh face normal
directions. The demo shows visualisations for a triangular, quadrilateral, and 
a pentagonal mesh. 
=#

fig = Figure(size=(1600,800))
cAlpha = RGBA(1.0, 1.0, 1.0,0.25)
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
        F,V = dodecahedron()
        invert_faces!(F)
        titleString="inverted pentagons"        
    end

    # Compute mesh face normals
    N = facenormal(F,V)    
    VN = simplexcenter(F,V)
    
    ax1=Axis3(fig[1, q], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString)
    hp1=poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=2,shading=FastShading,color=cAlpha, transparency=true, overdraw=false)
    # hpa=arrows!(ax1,VN,N,color=:blue)  
    hpa = dirplot(ax1,VN,N; color=:blue,linewidth=3,scaleval=1.0,style=:from)
end

fig