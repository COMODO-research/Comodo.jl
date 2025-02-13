using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

testCase = 1
if testCase == 1
    # Example geometry
    F,V = geosphere(3,1.0)
    VC = simplexcenter(F,V)
    F = [F[i] for i in findall(map(v-> v[3]>0,VC))] # Remove some faces
elseif testCase == 2
    F,V = geosphere(3,1.0)
    Vr = 2.0*randn(Point{3,Float64},5)
    V = [Vr;V]
    F = [f.+length(Vr) for f in F]
elseif testCase == 3
    F,V = geosphere(2,1.0)
    F = [F[i] for i in 1:2:length(F)]
end

indCheck = intersect(findall([v[1]>0 for v in V]),reduce(vcat,F))
Vn = deepcopy(V)

Fc,Vc,indFix = remove_unused_vertices(F,V)
indCheck_c = indFix[indCheck]

# Visualization
markersize = 10
markersize2 = 20
fig = Figure(size=(1200,1200))
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Mesh with unused vertices")
hp1 = poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=1,color=:white, strokecolor=:blue, shading = FastShading, transparency=false)
scatter!(Vn,markersize=markersize,color=:red)
scatter!(Vn[indCheck],markersize=markersize2,color=:blue)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Mesh with unused vertices removed")
hp2 = poly!(ax2,GeometryBasics.Mesh(Vc,Fc), strokewidth=1,color=:white, strokecolor=:blue, shading = FastShading, transparency=false)
scatter!(Vc,markersize=markersize,color=:red)
scatter!(Vc[indCheck_c],markersize=markersize2,color=:blue)

fig
