using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra

GLMakie.closeall()

nSteps = 50
xr,yr,zr = ntuple(_->range(0,4*pi,nSteps),3)

A = [triplyperiodicminimal((x,y,z), :G) for x in xr, y in yr, z in zr]
level = 0.0
cap = true

F1, V1 = getisosurface(A; x = xr, y = yr, z = zr, level = level, cap = cap, padValue=1e8)      

# Visualization
fig = Figure(size=(800,800))

ax1 = AxisGeom(fig[1, 1]; title="level=$level", limits=(minimum(xr),maximum(xr),minimum(yr),maximum(yr),minimum(zr),maximum(zr)))
hp1 = meshplot!(ax1, F1, V1, color=:brown, strokewidth=0.0)

stepRange1 = range(-1.5, 1.5, 50)
hSlider1 = Slider(fig[2, :], range = stepRange1, startvalue = level, linewidth=30)

function update_iso(A, level)
    F1, V1 = getisosurface(A; x=collect(xr), y=collect(yr), z=collect(zr), level=level, cap=cap, padValue=1e8)    
    F1, V1 = separate_vertices(F1, V1)    
    if !isempty(F1)
        hp1[1] = GeometryBasics.Mesh(V1, F1)
        hp1.visible=true
    else
        hp1.visible=false       
    end
end

on(hSlider1.value) do level
    A = [triplyperiodicminimal((x,y,z), :G) for x in xr, y in yr, z in zr]
    update_iso(A, level)
    ax1.title = "s=$s, level=$level"
end

screen = display(GLMakie.Screen(), fig)