using Comodo
using GLMakie
using GeometryBasics
using Printf

r = sqrt(3.0)
h = 0.5
F,V = dodecahedron(r;h=h)

## Visualize mesh
markersize = 25
strokewidth = 2 
strokecolor = :black

fig = Figure(size = (2000,1200))
ax1 = Axis3(fig[1:3, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title =  @sprintf "Pyritohedron h=%0.6f" h)

Fn,Vn = separate_vertices(F,V)

hp1 = poly!(ax1, GeometryBasics.Mesh(Vn,Fn), color=:white,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)
hp2 = scatter!(ax1, V,markersize=markersize,color=:red)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title =  "Cube")
Fs,Vs = cube(1.0);Fs,Vs = separate_vertices(Fs,Vs)
poly!(ax2, GeometryBasics.Mesh(Vs,Fs), color=:white,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)

ax3 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title =  "Dodecahedron")
Fs,Vs = dodecahedron(1.0);Fs,Vs = separate_vertices(Fs,Vs)
poly!(ax3, GeometryBasics.Mesh(Vs,Fs), color=:white,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)

ax4 = Axis3(fig[3, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title =  "Rhombic dodecahedron")
Fs,Vs = rhombicdodecahedron(1.0);Fs,Vs = separate_vertices(Fs,Vs)
poly!(ax4, GeometryBasics.Mesh(Vs,Fs), color=:white,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)

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