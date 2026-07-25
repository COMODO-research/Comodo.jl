using Comodo
using Comodo.GLMakie

r = 2.0 # Cylinder radius
h = 5.0 # Cylinder height

# Using splitting for regular cap meshes
n_set = (0, 1, 2)

# Loop over options and visualise
GLMakie.closeall()
cmap = Makie.Categorical(:Spectral) 

fig = Figure(size=(1200,800))

for (i,n) in enumerate(n_set)
    F, V, C = quadcylinder(r, h, n)

    # Visualise
    Fs, Vs = separate_vertices(F, V)
    Cs = simplex2vertexdata(Fs, C, Vs)
    ax1 = AxisGeom(fig[1, i], title = "n="*string(n))
    hp1 = meshplot!(ax1, Fs, Vs; color=Cs, colormap=cmap)
    # normalplot(ax1, F, V)
end

fig