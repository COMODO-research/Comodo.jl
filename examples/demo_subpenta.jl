using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of `subpenta` to refine a pentahedral mesh through 
splitting. 
=#

# Create example input mesh
r = 1.0 # Radius
n = 0
t = 1.0
nSteps = 2
direction=:positive
Fs,Vs = tridisc(r,n)
E, V = extrudefaces(Fs,Vs; extent=t, direction=direction, num_steps=nSteps)
F = element2faces(E)

# Define input parameters
nRefine = 1 # Number of refinement steps

# Create three variations i.e. direction 0, 1 and 2. 
Es0,Vs0 = subpenta(E, V, nRefine; direction=0)
Fs0 = element2faces(Es0)

Es1,Vs1 = subpenta(E, V, nRefine; direction=1)
Fs1 = element2faces(Es1)

Es2,Vs2 = subpenta(E, V, nRefine; direction=2)
Fs2 = element2faces(Es2)

# Visualisation
strokewidth = 2
linewidth = 3

fig = Figure(size=(1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Original mesh")
hp1 = meshplot!(ax1, F[1], V, strokewidth=strokewidth)
hp1 = meshplot!(ax1, F[2], V, strokewidth=strokewidth)

hpa = normalplot(ax1, F[1], V; color=:blue,linewidth=3)
hpa = normalplot(ax1, F[2], V; color=:blue,linewidth=3)

ax2 = AxisGeom(fig[1, 2], title = "Refined mesh, direction=0 (all)")
hp1 = meshplot!(ax2, Fs0[1], Vs0, strokewidth=strokewidth)
hp1 = meshplot!(ax2, Fs0[2], Vs0, strokewidth=strokewidth, transparency=false)

hpa = normalplot(ax2, Fs0[1], Vs0; color=:blue,linewidth=3)
hpa = normalplot(ax2, Fs0[2], Vs0; color=:blue,linewidth=3)

ax3 = AxisGeom(fig[2, 1], title = "Refined mesh, direction=1 (in plane)")
hp1 = meshplot!(ax3, Fs1[1], Vs1, strokewidth=strokewidth)
hp1 = meshplot!(ax3, Fs1[2], Vs1, strokewidth=strokewidth, transparency=false)

hpa = normalplot(ax3, Fs1[1], Vs1; color=:blue,linewidth=3)
hpa = normalplot(ax3, Fs1[2], Vs1; color=:blue,linewidth=3)

ax4 = AxisGeom(fig[2, 2], title = "Refined mesh, direction=2 (allong length)")
hp1 = meshplot!(ax4, Fs2[1], Vs2, strokewidth=strokewidth)
hp1 = meshplot!(ax4, Fs2[2], Vs2, strokewidth=strokewidth, transparency=false)

hpa = normalplot(ax4, Fs2[1], Vs2; color=:blue,linewidth=3)
hpa = normalplot(ax4, Fs2[2], Vs2; color=:blue,linewidth=3)

fig