using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Statistics

r = 10.0 # Radius
n = 3 # Number of refinement steps from cube
F, V, C = subquadsphere(n, r)

## Visualization
GLMakie.closeall()

Fs, Vs = separate_vertices(F, V)
Cs = simplex2vertexdata(Fs, C, Vs)

fig = Figure(size=(800,800))

ax1 = AxisGeom(fig[1, 1], title = "quadrangulated sphere", limits=(-r, r, -r, r, -r, r))
hp1 = meshplot!(ax1, Fs, Vs, strokewidth=strokewidth, color=Cs, colormap=cmap)
Colorbar(fig[1, 2], hp1)

stepRange = 0:6
hSlider = Slider(fig[2, :], range = stepRange, startvalue = 3, linewidth=30)

options1 = zip(["linear", "Catmull_Clark"], [:linear, :Catmull_Clark])
menu1 = Menu(fig[3,:], options = options1, default = "linear")

options2 = zip(["cube", "rhombicdodecahedron"], [:cube, :rhombicdodecahedron])
menu2 = Menu(fig[4,:], options = options2, default = "cube")

function updatePlot_subquadsphere(hSlider, menu1, menu2)
    n = hSlider.value[] 
    m = menu1.selection.val
    t = menu2.selection.val

    F, V, C = subquadsphere(n, r; method=m, template=t)
    Fs, Vs = separate_vertices(F,V)
    Cs = simplex2vertexdata(Fs,C)

    hp1[1] = GeometryBasics.Mesh(Vs, Fs)
    hp1.color = Cs
    ax1.title = "quadrangulated sphere, n=$n, method=$m, template=$t"
end

on(hSlider.value) do n
    updatePlot_subquadsphere(hSlider, menu1, menu2)
end

on(menu1.selection) do m
    updatePlot_subquadsphere(hSlider, menu1, menu2)
end

on(menu2.selection) do m
    updatePlot_subquadsphere(hSlider, menu1, menu2)
end

fig