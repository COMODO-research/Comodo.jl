using Comodo
using Comodo.GLMakie
using Comodo.Distributions

r = 2.0
N = 1000
P1 = rand_oncircle(r,N)
P2 = rand_incircle(r,N)
P3 = rand_onsphere(r,N)
P4 = rand_insphere(r,N)

# -----------------------------------------------------------------------------
# Visualisation

GLMakie.closeall()

F, V = geosphere(4, r)
V1 = circlepoints(r, 250)

fig = Figure(size=(1000, 1000))

ax1 = AxisGeom(fig[1, 1]; title="random points on circle", azimuth=-pi/2, elevation=pi/2)
lines!(ax1, V1, linewidth=6, color=:red)
hp1 = scatter!(ax1, P1, markersize=5, color=:black)

ax2 = AxisGeom(fig[1, 2]; title="random points in circle", azimuth=-pi/2, elevation=pi/2)
lines!(ax2, V1, linewidth=4, color=:red)
hp2 = scatter!(ax2, P2, markersize=5, color=:black)

ax3 = AxisGeom(fig[2, 1]; title="random points on sphere")
meshplot!(ax3, F, V; color=:red, strokewidth=0.0)
hp3 = scatter!(ax3, P3, markersize=5, color=:black, depth_shift=-0.01f0)

ax4 = AxisGeom(fig[2, 2]; title="random points in sphere")
hp4 = scatter!(ax4, P4, markersize=5, color=:black)
meshplot!(ax4, F, V; color=(:red, 0.2), transparency=true, strokewidth=0.0)

fig