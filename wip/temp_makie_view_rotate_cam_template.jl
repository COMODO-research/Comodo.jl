using Gibbon, GLMakie, GeometryBasics

r = 0.5 #radius
n = 3

M = platonicsolid(4,r)

fig = Figure()
scene = LScene(fig[1, 1])
cc = Makie.Camera3D(scene.scene, projectiontype = Makie.Perspective)

mesh!(scene,M)

center!(scene.scene)

fig