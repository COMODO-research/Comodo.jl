using Comodo
using GLMakie
using GeometryBasics
using Printf

n = 10
r = 2.0
F,V = ntrapezohedron(n,r)

## Visualize mesh
markersize = 25
strokewidth = 2 
strokecolor = :black

Fn,Vn = separate_vertices(F,V)

fig = Figure(size = (1200,1200))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "n-trapezohedron")
# ax1 = LScene(fig[1,1]); cc = Makie.Camera3D(ax1.scene, projectiontype = Makie.Orthographic)

hp1 = poly!(ax1, GeometryBasics.Mesh(Vn,Fn), color=:white,transparency=false,strokewidth=strokewidth,strokecolor=strokecolor,shading = FastShading)
hp2 = scatter!(ax1, V,markersize=markersize,color=:red)
# hp3 = normalplot(ax1,F,V; color = :green)

stepRange = 3:10
hSlider = Slider(fig[2,1], range = stepRange, startvalue = 0,linewidth=30)

on(hSlider.value) do n 
    F,V = ntrapezohedron(n,r)    
    Fn,Vn = separate_vertices(F,V)
    M = GeometryBasics.Mesh(Vn,Fn)
    hp1[1] = M
    hp2[1] = V    
    ax1.title =  @sprintf "n-trapezohedron n=%i" n  
end
slidercontrol(hSlider,ax1)

# fileName = comododir()*"/assets/temp/ntrapezohedron.mp4"
# slider2anim(fig,hSlider,fileName; backforth=true, duration=4)

fig
