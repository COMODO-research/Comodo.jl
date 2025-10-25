using Comodo
using Comodo.GLMakie
using Comodo.LinearAlgebra
using Comodo.Statistics
using Comodo.GeometryBasics
using Comodo.Rotations
using FileIO


testCase = 2

if testCase == 1
    F, V = geosphere(2, 1.0)
    V_mean = Point{3,Float64}(5.0, 2.0, 3.0)
    V = [Point{3,Float64}(4.0*v[1], 2.0*v[2], 0.75*v[3]) for v in V]

    # Define a rotation tensor using Euler angles
    Q = RotXYZ(0.25*π, -0.25*π, 0.25*pi)

    # Rotate the coordinates
    V = [Point{3, Float64}(Q*v) for v ∈ V] 
    V .+= V_mean
elseif testCase == 2
    # Loading a mesh
    fileName_mesh = joinpath(comododir(),"assets","stl","femur.stl")
    M = load(fileName_mesh)

    # Obtain mesh faces and vertices
    F = [TriangleFace{Int}(f) for f in faces(M)]
    V = coordinates(M)
    V = [v*1000 for v in V]
    F, V = mergevertices(F, V)
    V_mean = mean(V)
end 

S = pointsvd(V)

R = S.V # Orientations

Vr = [Point{3, Float64}(S.Vt*v) for v ∈ V] 

#Visualize mesh
GLMakie.closeall()

fig = Figure(size = (1600,800))

ax1 = AxisGeom(fig[1, 1], title = "SVD derived orientations")
hp1 = meshplot!(ax1, F, V; color=(:white, 0.5), strokewidth=0.0, transparency=true)
arrows3d!(ax1, V_mean, Point{3,Float64}(R[:,1]), color = :red, lengthscale = 100)
arrows3d!(ax1, V_mean, Point{3,Float64}(R[:,2]), color = :green, lengthscale = 100)
arrows3d!(ax1, V_mean, Point{3,Float64}(R[:,3]), color = :blue, lengthscale = 100)

ax1 = AxisGeom(fig[1, 2], title = "Rotated")
hp1 = meshplot!(ax1, F, Vr; strokewidth=0.0)

fig