using Comodo
using GLMakie
using GeometryBasics

#=
This demo shows the use of `lerp` for linear data interpolation. 
=#


# 1D curve interpolation 
x = range(0,9,9) # Data sites
y = 5.0*cos.(x.^2 ./ 9.0) # Data values
xi = range(-0.5,9.5,50) # Interpolation sites
yi = lerp(x,y,xi) # Linearly interpolate data 


# 3D curve interpolation based on 1D parameterisation
np = 10
t = range(0.0,2.0*π,np) # Parameterisation metric
V = [GeometryBasics.Point{3, Float64}(cos(t[i]),sin(t[i]),t[i]/(2.0*π)) for i ∈ eachindex(t)] # ND data, here 3D points
np_i = np*3
ti = range(minimum(t)-0.5,maximum(t)+0.5,np_i)
Vi = lerp(t,V,ti) # Linearly interpolate data 


# Visualization
fig = Figure(size = (1200,800))

ax1 = Axis(fig[1, 1], aspect = DataAspect(),title = "1D lerp interpolation")
hp1 = scatter!(ax1, x,y,markersize=25,color=:black)
scatter!(ax1, xi,yi,markersize=15,color=:red)
hp2 = lines!(ax1, xi,yi,linewidth=3,color=:red)
Legend(fig[1, 2],[hp1,hp2],["Raw","Interpolated"])

ax2 = Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "ND lerp interpolation")

hp1 = scatter!(ax2,V,markersize=25,color=:black)
scatter!(ax2, Vi,markersize=15,color=:red)
hp2 = lines!(ax2, Vi,linewidth=3,color=:red)
Legend(fig[1, 4],[hp1,hp2],["Raw","Interpolated"])

fig
