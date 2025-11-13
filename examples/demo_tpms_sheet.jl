using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra

GLMakie.closeall()

s = 0.0
level = 0.3
x = range(0.0, 4.0*pi, 50)
type=:G
F1, V1 = tpms_sheet(type, s; x=x, level=level, cap = true, padValue=1e8, side=:positive)
F2, V2 = tpms_sheet(type, s; x=x, level=level, cap = true, padValue=1e8, side=:negative)

# Visualization
fig = Figure(size=(1200,800))

ax1 = AxisGeom(fig[1, 1]; title="positive")
hp1 = meshplot!(ax1, F1, V1, color=:brown, strokewidth=0.0)

ax2 = AxisGeom(fig[1, 2]; title="negative")
hp2 = meshplot!(ax2, F2, V2, color=:white, strokewidth=0.0)

ax3 = AxisGeom(fig[1, 3]; title="s=$s, level=$level")
hp3 = meshplot!(ax3, F1, V1, color=:brown, strokewidth=0.0)
hp4 = meshplot!(ax3, F2, V2, color=(:white, 0.5), strokewidth=0.0, transparency=true)

stepRange1 = range(-1.5, 1.5, 50)
hSlider1 = Slider(fig[2, :], range = stepRange1, startvalue = level, linewidth=30)

stepRange2 = range(0.0, 1.5, 50)
hSlider2 = Slider(fig[3, :], range = stepRange2, startvalue = level, linewidth=30)

function update_iso_sheet(s, level)
    F1, V1 = tpms_sheet(type, s; x=x, level=level, cap = true, padValue=1e8, side=:positive)  
    F2, V2 = tpms_sheet(type, s; x=x, level=level, cap = true, padValue=1e8, side=:negative)  

    if !isempty(F1)
        hp1[1] = GeometryBasics.Mesh(V1, F1)
        hp1.visible=true
        hp2[1] = GeometryBasics.Mesh(V2, F2)
        hp2.visible=true
        hp3[1] = GeometryBasics.Mesh(V1, F1)
        hp3.visible=true
        hp4[1] = GeometryBasics.Mesh(V2, F2)
        hp4.visible=true
    else
        hp1.visible=false
        hp2.visible=false      
        hp3.visible=false
        hp4.visible=false   
    end
end

on(hSlider1.value) do s
    level = hSlider2.value[]
    update_iso_sheet(s, level)
    ax3.title = "s=$s, level=$level"
end

on(hSlider2.value) do level
    s = hSlider1.value[]
    update_iso_sheet(s, level)
    ax3.title = "s=$s, level=$level"
end

screen = display(GLMakie.Screen(), fig)