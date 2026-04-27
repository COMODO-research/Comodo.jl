# Quadrilateral Plate Mesh

This example demonstrates how to generate and visualize a quadrilateral plate mesh using Comodo.jl.

## Load packages

```julia
using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
```

## Generate plate mesh

```julia
plateDim = [20.0, 24.0]
plateElem = [11, 16]
orientation = :up

F, V, Eb, Cb = quadplate(plateDim, plateElem; orientation = orientation)
```

## Visualize mesh

```julia
GLMakie.closeall()

cmap = Makie.Categorical(:Spectral)

Ebs, Vbs = separate_vertices(Eb, V)
Cbs_V = simplex2vertexdata(Ebs, Cb)

fig = Figure(size = (1200, 800))

ax1 = AxisGeom( fig[1, 1],   title = "Quadrilateral mesh plate",   azimuth = -pi/2,  elevation = pi/2)

meshplot!(ax1, F, V)

edgeplot!(ax1,Ebs,Vbs;color = Cbs_V,linewidth = 6,colormap = cmap )

Colorbar(fig[1, 2])

fig
```

