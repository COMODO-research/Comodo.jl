using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.GLMakie.Colors
using Comodo.Statistics

#=
This demo shows the use of `hexahedronElement` to return a single hexahedral 
element. 
=#

w = 2.0 # Element width
e, V = hexahedronElement(w)

# Visualization
F = element2faces(e)
Fs, Vs = separate_vertices(F, V)

fig = Figure(size=(1700, 1000))
ax1 = AxisGeom(fig[1,1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "A hexahedral element")
hp1 = meshplot!(ax1, Fs, Vs, color=(:white, 0.25), strokewidth=3, strokecolor=:blue, transparency=true)
screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "testCase = $testCase")