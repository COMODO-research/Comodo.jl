using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.Statistics
using FileIO

#=
This demo shows the use of the `subtri_centre` function to sub-devide the faces
of a triangulated mesh iterative, by introducing a new central face point for 
each triangle, and by forming 3 new triangles for each triangle using the 
central point and the three triangle edges. 
=#

# Loading a mesh
fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
M = load(fileName_mesh)

# Obtain mesh faces and vertices
F = [TriangleFace{Int}(f) for f in faces(M)]
V = coordinates(M)
F, V = mergevertices(F, V) 

n = 1 # Number of splitting iterations
Fn, Vn = subtri_centre(F, V, n) # subdevide triangular mesh

# Visualisation
GLMakie.closeall()

fig = Figure(size = (1200,800))
ax1 = AxisGeom(fig[1, 1])
hp1 = meshplot!(ax1, F,  V, color=:white)

ax2 = AxisGeom(fig[1, 2], title="n = $n refinement steps")
hp2 = meshplot!(ax2, Fn,  Vn, color=:white, strokewidth=1.0)

stepRange = 0:1:5
hSlider = Slider(fig[2, :], range = stepRange, startvalue = 1,linewidth=30)

on(hSlider.value) do n    
    Fn,Vn = subtri_centre(F, V, n)     
    ax2.title = "n = $n refinement steps"
    hp2[1] = GeometryBasics.Mesh(Vn, Fn)
end

slidercontrol(hSlider, ax2)
fig
