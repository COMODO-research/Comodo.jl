using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Statistics
using FileIO

#=
This demo shows the use of `subtri_dual` to refine triangulated meshes. The
dual subtriangulation method is sometimes referred to as √3-subdevision (see 
also [1]). Following 2k subdevision steps each original triangle now represents 
9ᵏ triangles. The implementation here largely follows reference [1] however 
smoothing can be turned on/off using the `smooth` option (default is `true`). 
In addition, there is the option to contrain the boundary during smoothing by 
setting `constrain_boundary=true` (default `false`). Finally there is the option 
to split the boundary edges by setting `split_boundary=true` (default is 
`false`).      
# References    
[1] Kobbelt, √3-subdivision, 2000, SIGGRAPH 2000. 
=#

GLMakie.closeall()
r = 160.0 # Radius
n1 = 3
F,V = tridisc(r,n1)
split_boundary = true
smooth = true        


function subtri_dual_local!(F,V,B)
    F_tri, V_tri = subtri_dual(F[B], V, 1; smooth=smooth, split_boundary=false)
    F = F[.!B]
    append!(F, F_tri)
    # return mergevertices(F,V_tri)
    return F, V_tri
end

d = 100.0
B = [norm(mean(V[f]))<d for f in F]
F_tri, V_tri = subtri_dual_local!(F,V,B)

d = 50.0
B = [norm(mean(V_tri[f]))<d for f in F_tri]
F_tri, V_tri= subtri_dual_local!(F_tri, V_tri,B)

# n = 1
# F_tri, V_tri = subtri_dual(F[B], V, n; smooth=true, split_boundary=false)


## Visualization
strokewidth1 = 0.5
lineWidth = 4

Fs, Vs = separate_vertices(F,V)

fig = Figure(size=(1600,800))

ax1 = AxisGeom(fig[1, 1], title = "Original")
meshplot!(ax1, Fs, Vs, strokewidth=strokewidth1, strokecolor=:red, color=:white, transparency=false)

ax2 = AxisGeom(fig[1, 2], title = "Refined")
hp = meshplot!(ax2, F_tri, V_tri, strokewidth=strokewidth1, strokecolor=:red, color=:white, transparency=false)

screen = display(GLMakie.Screen(), fig)
