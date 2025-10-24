using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of `lerp` for linear data interpolation. 
=#

GLMakie.closeall()

for testCase = 1:4
    if testCase == 1
        # 1D curve interpolation 
        x = range(0, 9, 9) # Data sites
        y = 5.0*cos.(x.^2 ./ 5.0) # Data values
        xi = range(-1.5, 10.5, 100) # Interpolation sites                

        yi = lerp(x, y, xi) # Linearly interpolate data        
    elseif testCase == 2
        # 1D curve interpolation reversed curve 
        x = range(0, -9, 9) # Data sites
        y = 5.0*cos.(-x.^2 ./ 5.0) # Data values
        xi = range(1.5, -10.5, 100) # Interpolation sites                

        yi = lerp(x, y, xi) # Linearly interpolate data
    elseif testCase == 3 
        # 1D curve interpolation reversed curve 
        x = collect(range(0, -9, 9)) # Data sites
        y = 5.0*cos.(-x.^2 ./ 5.0) # Data values
        xi = collect(range(1.5, -10.5, 100)) # Interpolation sites                
        reverse!(x)
        reverse!(y)        
        yi = lerp(x, y, xi) # Linearly interpolate data

    elseif testCase == 4
        # 3D curve interpolation based on 1D parameterisation
        np = 10
        t = range(0.0,2.0*π,np) # Parameterisation metric
        V = [GeometryBasics.Point{3, Float64}(cos(t[i]),sin(t[i]),t[i]/(2.0*π)) for i ∈ eachindex(t)] # ND data, here 3D points
        np_i = np*3
        ti = range(minimum(t)-0.5, maximum(t)+0.5, np_i)        

        Vi = lerp(t, V, ti) # Linearly interpolate data
    end

    # Visualization
    fig = Figure(size = (1200,800))
    if testCase == 4
        ax1 = AxisGeom(fig[1, 1], title = "ND lerp interpolation")
        hp1 = scatter!(ax1, V, markersize=25, color=:black)
        hp2 = lines!(ax1, Vi, linewidth=3, color=:red)
        scatter!(ax1, Vi, markersize=15, color=:red)
    else
        ax1 = Axis(fig[1, 1], title = "ND lerp interpolation")
        hp1 = scatter!(ax1, x, y, markersize=25, color=:black)
        hp2 = lines!(ax1, xi, yi, linewidth=3, color=:red)
        scatter!(ax1, xi, yi, markersize=15, color=:red)
    end
    Legend(fig[1, 2],[hp1,hp2],["Raw","Interpolated"])

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")

end