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