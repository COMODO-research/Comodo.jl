using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

r = 1.0 #radius
F,V = platonicsolid(4,r)

n = 1
Fn,Vn = subtri(F,V,n; method=:Loop) 

# Visualisation
GLMakie.closeall()

fig = Figure(size=(800,800))
ax = AxisGeom(fig[1, 1], title = "n = " * string(n) * " refinement steps")
hp1 = edgeplot!(ax, F, V, linewidth=3, color=:red)
hp2 = meshplot!(ax, Fn, Vn)

stepRange = 0:1:5
hSlider = Slider(fig[2, 1], range = stepRange, startvalue = 1,linewidth=30)

on(hSlider.value) do n    
    Fn,Vn = subtri(F,V,n; method=:Loop)     
    ax.title = "n = " * string(n) * " refinement steps"
    hp2[1] = GeometryBasics.Mesh(Vn, Fn)
end

fig