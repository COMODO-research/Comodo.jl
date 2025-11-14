using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

GLMakie.closeall()

for testCase = 1:3
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

    Fc = deepcopy(F)
    Vc = deepcopy(V)
    indFix = remove_unused_vertices!(Fc,Vc)
    indCheck_c = indFix[indCheck]

    # Visualization
    markersize = 10
    markersize2 = 20
    fig = Figure(size=(1200,1200))
    ax1 = AxisGeom(fig[1, 1], title = "Mesh with unused vertices")
    hp1 = meshplot!(ax1, F, V)
    scatter!(Vn, markersize=markersize, color=:red)
    scatter!(Vn[indCheck], markersize=markersize2, color=:blue)

    ax2 = AxisGeom(fig[1, 2], title = "Mesh with unused vertices removed")
    hp2 = meshplot!(ax2, Fc , Vc)
    scatter!(Vc, markersize=markersize, color=:red)
    scatter!(Vc[indCheck_c], markersize=markersize2, color=:blue)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end