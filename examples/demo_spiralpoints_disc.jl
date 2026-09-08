using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

N = 75
V = spiralpoints_disc(N)

# Visualisation
GLMakie.closeall()

r = 1.0 # radius
n1 = 5
F1,V1 = tridisc(r, n1) # Create disc mesh 
N1 = vertexnormal(F1, V1)

function morphDisc(V1, N1, V)
    C1 = mindist(V1, V)
    dMax = maximum(C1)
    D1 = [r/3.0*(d/dMax)^2 for d in C1]    
    V1m = [v-d*n for (v, n, d) in zip(V1, N1, D1)]
    return V1m, C1
end

V1m, D1 = morphDisc(V1, N1, V)

cmap = :viridis # Makie.Reverse(:Spectral) #(:viridis, 0.75)
fig = Figure(size = (1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Spiral points")
hp1 = meshplot!(ax1, F1, V1m, color=D1, strokewidth=0.0, transparency=false, colormap=cmap)
hp2 = scatter!(ax1, V, color=:black, markersize=10, depth_shift=-0.01f0)

toggle = Toggle(fig[2, :][1, 1], active = false, length=64, markersize=32)

stepRange1 = 1:2*N
hSlider1 = Slider(fig[2, :][1, 2], range = stepRange1, startvalue = N,linewidth=32)

on(hSlider1.value) do N
    V = spiralpoints_disc(N)
    V1m, D1 = morphDisc(V1, N1, V)   
    hp2[1] = V
    hp1[1] = GeometryBasics.Mesh(V1m, F1)
    hp1.color = D1
end

on(fig.scene.events.tick) do tick 
    toggle.active[] || return     
    sliderRange = hSlider1.range[] # Get slider range
    rangeLength = length(sliderRange) # Number of possible steps 
    sliderIndex = hSlider1.selected_index[]

    if sliderIndex < rangeLength                       
        sliderIndex += 1                                 
    else
        sliderIndex = 1
    end 
    hSlider1.selected_index = sliderIndex         
end

fig