
using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
 
GLMakie.closeall()

h = 80.0
n = 0
r = 25.0
nh = 0

E, V, F, Fb, Cb = hexcylinder(r, h, n; nh=nh, direction=:both)  

Fbs, Vbs = separate_vertices(Fb,V)
Cbs = simplex2vertexdata(Fbs,Cb)

# Visualisation
cmap = Makie.Categorical(:Spectral) 
fig = Figure(size=(1200,1000))
ax1 = AxisGeom(fig[1, 1], title="Hexahedral mesh of a cylinder, n=$n refinement steps")
hp1 = meshplot!(ax1,Fbs,Vbs; color=Cbs, colormap=cmap, strokecolor=:white, strokewidth=1)
Colorbar(fig[1, 1][1, 2], hp1)

stepRange = 0:4
hSlider1 = Slider(fig[2, :], range = stepRange, startvalue = 0, linewidth=30)

on(hSlider1.value) do n     
    h = hSlider2.value[]
    r = hSlider3.value[]
    E, V, F, Fb, Cb = hexcylinder(r, h, n; nh=nh)  
    Fbs, Vbs = separate_vertices(Fb,V)
    Cbs = simplex2vertexdata(Fbs,Cb)

    hp1[1] = GeometryBasics.Mesh(Vbs,Fbs)    
    hp1.color = Cbs
    ax1.title = "Hexahedral mesh of a cylinder, n=$n refinement steps"
end

stepRange = range(0.01,h,20)
hSlider2 = Slider(fig[3, :], range = stepRange, startvalue = h, linewidth=30)

on(hSlider2.value) do h     
    n = hSlider1.value[]
    r = hSlider3.value[]
    E, V, F, Fb, Cb = hexcylinder(r, h, n; nh=nh)  
    Fbs, Vbs = separate_vertices(Fb,V)
    Cbs = simplex2vertexdata(Fbs,Cb)

    hp1[1] = GeometryBasics.Mesh(Vbs,Fbs)    
    hp1.color = Cbs
    ax1.title = "Hexahedral mesh of a cylinder, n=$n refinement steps"
end

stepRange = range(0.01,r,20)
hSlider3 = Slider(fig[4, :], range = stepRange, startvalue = r, linewidth=30)

on(hSlider3.value) do r     
    n = hSlider1.value[]
    h = hSlider2.value[]
    E, V, F, Fb, Cb = hexcylinder(r, h, n; nh=nh, direction=:both)
    Fbs, Vbs = separate_vertices(Fb,V)
    Cbs = simplex2vertexdata(Fbs,Cb)

    hp1[1] = GeometryBasics.Mesh(Vbs,Fbs)    
    hp1.color = Cbs
    ax1.title = "Hexahedral mesh of a cylinder, n=$n refinement steps"
end


screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "Hexahedral mesh of a cylinder")

# fileName = "/home/kevin/Videos/hexcylinder.mp4"
# slider2anim(fig,hSlider,fileName; backforth=true, duration=4)