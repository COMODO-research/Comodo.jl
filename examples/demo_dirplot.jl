using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.GLMakie.Colors

#=
This demo shows the use of the `dirplot` function to visualize directional data. 
=#

GLMakie.closeall()

cAlpha = RGBA(1.0, 1.0, 1.0,0.25)
for testCase = 1:2
    if testCase == 1 
        P = Point{3,Float64}(0.0,0.0,0.0)
        N = Point{3,Float64}(0.0,0.0,1.0)

        Pv = Point{3,Float64}(0.0,0.0,0.0)
        Nv = Vec{3,Float64}(0.0,0.0,1.0)

        # Visualisation
        fig = Figure(size=(1600,800))
        styleSet = (:to,:from,:through)
        for i in eachindex(styleSet)
            ax1 = Axis3(fig[1, i], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = string(styleSet[i]))        
            hp1 = scatter!(ax1,P,markersize=35,color=:red)
            hpa = dirplot(ax1,P,N; color=:blue,linewidth=3,scaleval=1.0,style=styleSet[i])

            ax2 = Axis3(fig[2, i], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = string(styleSet[i]))        
            hp1 = scatter!(ax2,Pv,markersize=35,color=:red)
            hpa = dirplot(ax2,Pv,Nv; color=:blue,linewidth=3,scaleval=1.0,style=styleSet[i])
        end
        screen = display(GLMakie.Screen(), fig)
        GLMakie.set_title!(screen, "testCase = $testCase")        
    elseif testCase == 2 
        F,V = icosahedron()
        N = vertexnormal(F,V; weighting=:size)

        # Visualisation
        fig = Figure(size=(1600,800))
        styleSet = (:to,:from,:through)
        for i in eachindex(styleSet)
            ax1 = Axis3(fig[1, i], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = string(styleSet[i]))
            hp1 = poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=1,shading=FastShading,color=cAlpha, transparency=true, overdraw=false)        
            hpa = dirplot(ax1,V,N; color=:blue,linewidth=3,scaleval=1.0,style=styleSet[i])
        end
        screen = display(GLMakie.Screen(), fig)
        GLMakie.set_title!(screen, "testCase = $testCase")
    end
end