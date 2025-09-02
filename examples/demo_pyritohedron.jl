using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Printf

r = sqrt(3.0)
h = 0.5
F,V = dodecahedron(r;h=h)

## Visualize mesh
GLMakie.closeall()

markersize = 25
strokewidth = 2 
strokecolor = :black

fig = Figure(size = (2000,1200))
ax1 = AxisGeom(fig[1:3, 1], title =  @sprintf "Pyritohedron h=%0.6f" h)

Fn,Vn = separate_vertices(F,V)

hp1 = meshplot!(ax1, Fn, Vn, strokewidth=strokewidth,strokecolor=strokecolor)
hp2 = scatter!(ax1, V,markersize=markersize,color=:red)

ax2 = AxisGeom(fig[1, 2], title =  "Cube")
Fs,Vs = cube(1.0);Fs,Vs = separate_vertices(Fs,Vs)
meshplot!(ax2, Fs, Vs, strokewidth=strokewidth,strokecolor=strokecolor)

ax3 = AxisGeom(fig[2, 2], title =  "Dodecahedron")
Fs,Vs = dodecahedron(1.0);Fs,Vs = separate_vertices(Fs,Vs)
meshplot!(ax3, Fs, Vs, strokewidth=strokewidth,strokecolor=strokecolor)

ax4 = AxisGeom(fig[3, 2], title =  "Rhombic dodecahedron")
Fs,Vs = rhombicdodecahedron(1.0);Fs,Vs = separate_vertices(Fs,Vs)
meshplot!(ax4, Fs, Vs, strokewidth=strokewidth,strokecolor=strokecolor)

stepRange = collect(range(0.00001,0.9999,100))
hSlider = Slider(fig[4, :], range = stepRange, startvalue = 0,linewidth=30)

on(hSlider.value) do h 
    F,V = dodecahedron(r;h=h)    
    Fn,Vn = separate_vertices(F,V)
    M = GeometryBasics.Mesh(Vn,Fn)
    hp1[1] = M
    hp2[1] = V    
    ax1.title =  @sprintf "Pyritohedron h=%0.6f" h  
end
slidercontrol(hSlider,ax1)

# fileName = comododir()*"/assets/temp/pyritohedron.mp4"
# slider2anim(fig,hSlider,fileName; backforth=true, duration=4)

fig