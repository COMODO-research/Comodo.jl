
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
        colorNow=:red
    elseif q==2
        M=cube()
        titleString="quadrilaterals"
        colorNow=:green
    elseif q==3
        M=dodecahedron()
        titleString="pentagons"
        colorNow=:blue
    end

    # Compute mesh normals
    N = facenormal(M)
    VN = simplexcenter(faces(M),coordinates(M))
    
    ax1=Axis3(fig[1, q], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString)
    hp1=poly!(ax1,M, strokewidth=3,shading=FastShading,color=colorNow, transparency=true, overdraw=false)
    hpa=arrows!(ax1,VN,N)
end

fig
