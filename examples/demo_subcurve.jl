using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of `lerp` for linear data interpolation. 
=#

# 1D curve interpolation 
t = range(0.0,2*π,10) 
V = [Point{3,Float64}(ti,cos(ti),0.0) for ti in t] # Data values

n1 = 2
Vn1 = subcurve(V,n1)

n2 = 5
Vn2 = subcurve(V,n2)

# Visualization
markersize1 = 25
markersize2 = 15
linewidth1 = 4
linewidth2 = 2

fig = Figure(size = (1200,800))

ax1 = Axis(fig[1, 1], aspect = DataAspect(),title = "curve subdivision with $n1 added points")
hp1 = scatter!(ax1,V,markersize=markersize1,color=:red)
hp2 = lines!(ax1, V,linewidth=linewidth1,color=:red)
hp3 = scatter!(ax1,Vn1,markersize=markersize2,color=:black)
hp4 = lines!(ax1, Vn1,linewidth=linewidth2,color=:black)
Legend(fig[1, 2],[hp2,hp4],["Raw","Refined"])

ax2 = Axis(fig[2, 1], aspect = DataAspect(),title = "curve subdivision with $n2 added points")
hp1 = scatter!(ax2,V,markersize=markersize1,color=:red)
hp2 = lines!(ax2, V,linewidth=linewidth1,color=:red)
hp3 = scatter!(ax2,Vn2,markersize=markersize2,color=:black)
hp4 = lines!(ax2, Vn2,linewidth=linewidth2,color=:black)
Legend(fig[2, 2],[hp2,hp4],["Raw","Refined"])

fig