using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Rotations
using Comodo.GLMakie.Colors

# Example curves
r = 1.0
nc = 16
t = range(2.0*π-(2.0*π/nc),0,nc)
Vc = [Point3{Float64}(2.0+cos(tt),0.0,sin(tt)) for tt ∈ t]
n = Vec{3, Float64}(0.0,0.0,1.0)
num_steps = 25
close_loop = true


F,V = revolvecurve(Vc; extent = (2*pi - (2*pi/num_steps)), direction = :negative, 
n = n, num_steps = num_steps, periodicity = (true,true), face_type = :quad)

VF = simplexcenter(F,V)
z = [v[3] for v in VF]
L1 = z .< 0 
Ft = F[.!L1]
indTop = elements2indices(Ft)
indRand = unique(rand(1:length(indTop),200))
Nt = vertexnormal(Ft,V)          

## Visualization
GLMakie.closeall()

markersize = 0.1
doughColor = RGB(1.0, 0.74, 0.5) 
icingColor = RGB(145/255, 85/255, 77/255) 

fig = Figure(size=(1200,1200))

ax1 = AxisGeom(fig[1, 1], title = "Donut")
hp1 = meshplot!(ax1, F[L1], V, color = doughColor, strokewidth=0.0)
hp2 = meshplot!(ax1, F[.!L1], V, color = icingColor, strokewidth=0.0)

for i in indRand
    Fr,Vr = geosphere(2,0.07)
    Q = rotation_between(Nt[indTop[i]],Vec3(0.0,0.0,1.0)) # The rotation between the current normal and the z-axis
    Q2 = RotXYZ(0.0,0.0,rand(1)[1]*2*pi)
    Vr = [Q* Q2* Point{3,Float64}(v[1]*3,v[2],v[3]) for v in Vr]
    Vr .+= V[indTop[i]]
    meshplot!(ax1, Fr, Vr, color=rand(RGB,1)[1], strokewidth=0.0)
end
fig
