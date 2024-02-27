using Comodo
using GeometryBasics
using GLMakie

"""
This demo shows the use of sweep_loft
"""

# Define guide curve
P = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
P[1 ] = GeometryBasics.Point{3, Float64}( 0.0, 0.0, 0.0)
P[2 ] = GeometryBasics.Point{3, Float64}( 1.0, 0.0, 0.0)
P[3 ] = GeometryBasics.Point{3, Float64}( 1.0, 1.0, 0.0)
P[4 ] = GeometryBasics.Point{3, Float64}( 1.0, 1.0, 1.0)

n = 25 # Number of points
Vc = nbezier(P,n) # Get Bezier fit points


np = 8
t = range(0.0,2.0*π-(2.0*π/np),np)

# Define start entities
p1 = Vc[1]
n1 = Vec3{Float64}(1.0,0.0,0.0)
r1 = 0.1
v1 = [GeometryBasics.Point{3, Float64}(p1[1],p1[2]+r1*cos(tt),p1[3]+r1*sin(tt)) for tt ∈ t]

# Define end entities
p2 = Vc[end]
n2 = Vec3{Float64}(0.0,0.0,1.0)
r2 = 0.2
v2 = [GeometryBasics.Point{3, Float64}(p2[1]+r2*cos(tt),p2[2]+r2*sin(tt),p2[3]) for tt ∈ t]


# Visualization
fig = Figure(size = (800,800))
ax = Axis3(fig[1, 1],aspect = :data)

# hp1 = scatter!(ax, P,markersize=25,color=:black)
hp21 = scatter!(ax, Vc,markersize=15,color=:black)
hp22 = lines!(ax, Vc,linewidth=3,color=:black)

scatter!(ax, [p1,p2] ,markersize=25,color=:blue)

hp31 = scatter!(ax, v1,markersize=15,color=:green)
hp32 = lines!(ax, v1,linewidth=3,color=:green)

hp41 = scatter!(ax, v2,markersize=15,color=:green)
hp42 = lines!(ax, v2,linewidth=3,color=:green)

fig
