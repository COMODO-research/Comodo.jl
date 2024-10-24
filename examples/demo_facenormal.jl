using Comodo
using GLMakie
using GeometryBasics

#=
This demo shows the use of the `meshnormal` function to obtain mesh face normal
directions. The demo shows visualisations for a triangular, quadrilateral, and 
a pentagonal mesh. 
=#

fig = Figure(size=(1600,800))

for q=1:1:3
    if q==1
        F,V = icosahedron()
        titleString="triangles"        
    elseif q==2
        F,V = cube()
        titleString="quadrilaterals"        
    elseif q==3
        F,V = dodecahedron()
        titleString="pentagons"        
    end

    # Compute mesh face normals
    N = facenormal(F,V)    
    VN = simplexcenter(F,V)
    
    ax1=Axis3(fig[1, q], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString)
    hp1=poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=3,shading=FastShading,color=:white, transparency=true, overdraw=false)
    # hpa=arrows!(ax1,VN,N,color=:blue)  
    hpa = dirplot(ax1,VN,N; color=:blue,linewidth=3,scaleval=1.0,style=:from)
end

fig
