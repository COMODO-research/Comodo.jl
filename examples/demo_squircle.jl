using Comodo
using GeometryBasics
using GLMakie
using Printf

τ = 0.5
n = 4*200 
r = 1.0
a_tol = 1e-6



V = squircle(r,n,τ; atol=a_tol)

# Visualization

markersize = 6
linewidth = 10
nSlider = 150
cMap = cgrad(:Spectral, nSlider, categorical = true)#Makie.Reverse(:Spectral)

fig = Figure(size = (1600,1600))
# ax = Axis3(fig[1, 1],aspect = :data,title="Squircle, τ="*string(τ))
ax = Axis(fig[1, 1],aspect = DataAspect(),title="Squircle, τ="*@sprintf("%.3f",τ),titlesize=50)

stepRange1 = range(0.0,1.0,nSlider)

for (i,stepVal) in enumerate(stepRange1)
    V = squircle(r,n,stepVal; atol=a_tol)
    lines!(ax, V,linewidth=0.25,color=cMap[i],transparency=true,depth_shift=0.01f0)
end

hSlider1 = Slider(fig[2, 1], range = stepRange1, startvalue = 0,linewidth=30)

hp1 = lines!(ax, V,linewidth=linewidth,color=:black,depth_shift=-0.01f0)

on(hSlider1.value) do stepVal        
    hp1[1] = squircle(r,n,stepVal; atol=a_tol)
    ax.title = "Squircle, τ="*@sprintf("%.3f",stepVal)
end

Colorbar(fig[1, 2], limits = (0, 1), colormap = cMap)

fig

# fileName = comododir()*"/assets/img/squircle.mp4"
# slider2anim(fig,hSlider1,fileName; backforth=true, duration=3)