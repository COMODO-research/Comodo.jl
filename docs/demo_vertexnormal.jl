
using Gibbon
using GLMakie

"""
This demo shows the use of the `meshnormal` function to obtain mesh face normal
directions. The demo shows visualisations for a triangular, quadrilateral, and 
a pentagonal mesh. 
"""

fig = Figure(size=(1600,800))

for q=1:1:3
    if q==1
        M=icosahedron()
        titleString="triangles"        
    elseif q==2
        M=cube()
        titleString="quadrilaterals"        
    elseif q==3
        M=dodecahedron()
        titleString="pentagons"        
    end

    # Compute mesh vertex normals
    NV = vertexnormal(M)
    V = coordinates(M)

    ax1=Axis3(fig[1, q], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString)
    hp1=poly!(ax1,M, strokewidth=3,shading=FastShading,color=:white, transparency=true, overdraw=false)    
    hpa=arrows!(ax1,V,NV,color=:red)
end

fig
