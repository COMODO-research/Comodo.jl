using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Statistics

r = 10.0/pi # Radius

n1 = 0 # Number of refinement steps from cube
Fn1,Vn1 = subquadsphere(n1,1)

n2 = 1 # Number of refinement steps from cube
Fn2,Vn2 = subquadsphere(n2,r)

n3 = 2 # Number of refinement steps from cube
Fn3,Vn3 = subquadsphere(n3,r)

n4 = 3 # Number of refinement steps from cube
Fn4,Vn4 = subquadsphere(n4,r)

pointSpacing = 0.25

## Visualization
GLMakie.closeall()

strokewidth = 2

fig = Figure(size=(800,800))

ax1 = AxisGeom(fig[1, 1], title = "Refined n=$n1")
meshplot!(ax1, Fn1, Vn1, strokewidth=strokewidth)

ax2 = AxisGeom(fig[1, 2], title = "Refined n=$n2")
meshplot!(ax2, Fn2, Vn2, strokewidth=strokewidth)

ax3 = AxisGeom(fig[2, 1], title = "Refined n=$n3")
meshplot!(ax3, Fn3, Vn3, strokewidth=strokewidth)

ax4 = AxisGeom(fig[2, 2], title = "Refined n=$n4")
meshplot!(ax4, Fn4, Vn4, strokewidth=strokewidth)

fig