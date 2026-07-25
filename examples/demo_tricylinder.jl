using Comodo
using Comodo.GLMakie

r = 2.0 # Cylinder radius
h = 8.0 # Cylinder height

# Using splitting for regular cap meshes
n_set = (0, 1, 2, 2)
pointSpacing_set = (1.0, 0.5, 0.3, 0.3)
face_type_set = (:tri, :tri, :tri, :forwardslash)

# Loop over options and visualise
GLMakie.closeall()
cmap = Makie.Categorical(:Spectral) 

fig = Figure(size=(1200,800))

for (i,n) in enumerate(n_set)
    F, V, C = tricylinder(r, h, n; face_type=face_type_set[i])

    # Visualise
    Fs, Vs = separate_vertices(F, V)
    Cs = simplex2vertexdata(Fs, C, Vs)
    ax1 = AxisGeom(fig[1, i], title = "n="*string(n)*", face_type=" * string(face_type_set[i]))
    hp1 = meshplot!(ax1, Fs, Vs; color=Cs, colormap=cmap)
    # normalplot(ax1, F, V)
end

for (i,pointSpacing) in enumerate(pointSpacing_set)
    F, V, C = tricylinder(r, h, pointSpacing; face_type=face_type_set[i])

    # Visualise
    Fs, Vs = separate_vertices(F, V)
    Cs = simplex2vertexdata(Fs, C, Vs)
    ax1 = AxisGeom(fig[2, i], title = "pointSpacing="*string(pointSpacing)*", face_type=" * string(face_type_set[i]))
    hp1 = meshplot!(ax1, Fs, Vs; color=Cs, colormap=cmap)
    # normalplot(ax1, F, V)
end

fig