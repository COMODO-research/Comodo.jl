using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of the `normalplot` function to visualise the normal 
directions for a surface mesh. 

    Full syntax example: 
    hpa = normalplot(ax1,F,V; type_flag=:face, color=:blue,linewidth=3,scaleval=1.0)
=#


type_flag_set = (:face,:face,:vertex)
fig = Figure(size=(1600,800))
GLMakie.closeall()

for q=1:1:3
    type_flag = type_flag_set[q]
    if q==1
        F,V = icosahedron()
        titleString="triangles, surface normals, type_flag=$type_flag"        
    elseif q==2
        F,V = cube()
        titleString="quadrilaterals, surface normals, type_flag=$type_flag"           
    elseif q==3
        F,V = dodecahedron()
        titleString="pentagons, surface normals, type_flag=$type_flag"            
    end
    
    ax1 = AxisGeom(fig[1, q], title = titleString)
    hp1 = meshplot!(ax1, F, V, strokewidth=3)
    # hpa=arrows!(ax1,VN,N,color=:blue)  
    # type_flag=:face, color=:black,linewidth=3,scaleval=nothing
    hpa = normalplot(ax1,F,V; type_flag=type_flag, color=:blue,linewidth=3,scaleval=1.0)
end

fig