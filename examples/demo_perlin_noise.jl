using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Random

Random.seed!(1) # To ensure repeatable demo behaviour

# Define grid
size_grid = (25,30) # grid size
sampleFactor = 30 # Pixel sample factor wrt grid

M = perlin_noise(size_grid, sampleFactor; type=:Perlin)

# Visualization
fig = Figure(size=(800,800))
size_image = size(M)
ax1 = Axis(fig[1, 1], aspect = DataAspect(), title = "Perlin noise",limits=(-sampleFactor,size_image[1]+sampleFactor,-sampleFactor,size_image[2]+sampleFactor) )
hm = image!(ax1, M, interpolate = false, colormap = Makie.Reverse(:Spectral), colorrange=(-0.5,0.5))
Colorbar(fig[1, 2], hm)
fig