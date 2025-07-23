using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
 
GLMakie.closeall()

n = 12
r = 25.0

F, V = pizza(r, n; dir=:acw)  

# Visualisation
fig = Figure(size=(1200,1000))
ax1 = AxisGeom(fig[1, 1], title="Pizza mesh, n=$n slices")
hp1 = meshplot!(ax1,F,V; color=:red, strokecolor=:white, strokewidth=1)

stepRange = 3:24
hSlider1 = Slider(fig[2, :], range = stepRange, startvalue = 12, linewidth=30)

on(hSlider1.value) do n     
    F, V = pizza(r, n; dir=:acw) 
    hp1[1] = GeometryBasics.Mesh(V,F)        
    ax1.title = "Pizza mesh, n=$n slices"
end

screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "Hexahedral mesh of a cylinder")