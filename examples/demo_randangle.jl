using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Random

Random.seed!(1) # To ensure repeatable demo behaviour

A1 = randangle(5)
println("Vector: ")
display(A1)


siz = (25,25)
A2 = randangle(siz)
println("2D array/matrix: ")
display(A2)

siz = (5,4,3)
A3 = randangle(siz)
println("3D array: ")
display(A3)

# # Visualization
# fig = Figure(size=(800,800))
# size_image = size(M)
# ax1 = Axis(fig[1, 1], aspect = DataAspect(), title = "An array of random angular data")
# hm = image!(ax1, A, interpolate = false, colormap = Makie.Reverse(:Spectral), colorrange=(-pi,pi))
# Colorbar(fig[1, 2], hm)
# fig