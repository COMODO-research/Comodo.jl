using GeometryBasics # For point and mesh format
using GLMakie
using Gibbon

r = 1.0 
ϕ = Base.MathConstants.golden # (1.0+sqrt(5.0))/2.0, Golden ratio
s = r/sqrt(3.0)
t = ϕ*s    
w = (ϕ-1.0)*s

r = 0.5 #radius
n = 3

M = platonicsolid(4,r)
V = coordinates(M)
F = faces(M)

Fn,Vn = subtri(F,V,n)


Dn = minDist(Vn,V; getIndex = false)

###
# Visualize mesh
fig = Figure()

ax=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")

hp=poly!(ax,GeometryBasics.Mesh(Vn,Fn),strokewidth=1,color=Dn, transparency=false, overdraw=false)
Colorbar(fig[1, 2], hp)
fig


# # Visualize 
# fig1 = Figure()

# ax = Axis3(fig1[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")

# scatter!(ax, V1,markersize=25,color=D)
# scatter!(ax, V2,markersize=25,color=:black)

# fig1