using Comodo
using Comodo.GLMakie
using Comodo.LinearAlgebra

#=
This demo shows the use of the gradient function to compute the gradient of an
array allong a certain direction. 
=#

GLMakie.closeall()

for testCase = 1:3
    if testCase == 1
        t = range(0.0,2*pi,12)
        f = sin.(t)
        df_dt = cos.(t)
        dims = 1
        Δx = t[2]-t[1]
        G = Vector{Float64}(undef,length(f))
        for i in eachindex(f)
            G[i] = gradient(f, Δx, i, dims)[1]
        end

        # Visualisation
        fig = Figure(size=(800,800))
        ax1 = Axis(fig[1, 1], title = "Input and gradient curve")
        hp1 = lines!(ax1, t, f, color=:blue, linewidth=2)
        hp2 = lines!(ax1, t, df_dt, color=:green, linewidth=2)
        hp3 = lines!(ax1, t, G, color=:red, linewidth=2)
        Legend(fig[1,2],[hp1, hp2, hp3],["f(t)", "df/dt", "gradient of f"])
        screen = display(GLMakie.Screen(), fig)
        GLMakie.set_title!(screen, "testCase = $testCase")
    elseif testCase == 2
        x = range(0.0, 2*pi, 25)
        y = range(0.0, 2*pi, 25)
        f = Matrix{Float64}(undef, length(x), length(y))
        df_dx = Matrix{Float64}(undef, length(x), length(y))
        df_dy = Matrix{Float64}(undef, length(x), length(y))
        for (i, x) in enumerate(x)
            for (j, y) in enumerate(y)
                f[i,j] = sin(x)+cos(y)
                df_dx[i, j] = cos(x) 
                df_dy[i, j] = -sin(y) 
            end
        end
        
        dims = [1,2]
        Δx = [x[2]-x[1], y[2]-y[1]]
        Gx = Matrix{Float64}(undef,size(f))
        Gy = Matrix{Float64}(undef,size(f))

        for i in eachindex(f)
            g = gradient(f, Δx, i, dims)
            Gx[i] = g[1]
            Gy[i] = g[2]
        end

        # Visualisation
        fig = Figure(size=(800,800))
        ax1 = AxisGeom(fig[1, 1], title = "Input f")
        hp1 = surface!(ax1, x, y, f)

        ax2 = AxisGeom(fig[1, 2][1, 1], title = "df/dx")
        hp2 = surface!(ax2, x, y, df_dx)

        ax3 = AxisGeom(fig[1, 2][2, 1], title = "df/dy")
        hp3 = surface!(ax3, x, y, df_dy)

        ax4 = AxisGeom(fig[1, 3][1, 1], title = "gradient x")
        hp4 = surface!(ax4, x, y, Gx)

        ax5 = AxisGeom(fig[1, 3][2, 1], title = "gradient y")
        hp5 = surface!(ax5, x, y, Gy)

        screen = display(GLMakie.Screen(), fig)
        GLMakie.set_title!(screen, "testCase = $testCase")

        # G = gradient(f, Δx, ind, dims)
    elseif testCase == 3
        # Create example 3D image data 
        nSteps = 25
        w = 4.77
        xr,yr,zr = ntuple(_->range(-w, w,nSteps),3)
        function thisImage(x,y,z)
            phi=(1+sqrt(5))/2
            return 1.0/6.0 *(2.0 - (cos(x + phi*y) + cos(x - phi*y) + cos(y + phi*z) + cos(y - phi*z) + cos(z - phi*x) + cos(z + phi*x)))        
        end
        I = [thisImage(x,y,z) for x in xr, y in yr, z in zr] 
        level = 0.0
        B = [i<=2.0 for i in I]
        indB = findall(B)

        P0 = Vector{Point{3, Float64}}(undef, length(indB))
        P = Vector{Point{3, Float64}}(undef, length(indB))
        C = Vector{Float64}(undef, length(indB))
        Δx = xr[2] - xr[1]
        Δy = yr[2] - yr[1]
        Δz = zr[2] - zr[1]
        for (j, i) in enumerate(indB)
            c = CartesianIndices(I)[i]
            x = gradient(I, Δx, i, 1)[1]
            y = gradient(I, Δy, i, 2)[1]
            z = gradient(I, Δz, i, 3)[1]
            p = Point{3, Float64}(x, y, z)
            P[j] = p
            P0[j] = Point{3, Float64}(xr[c[1]], yr[c[2]], zr[c[3]])
            C[j] = norm(p)
        end

        # Visualisation
        F, V = getisosurface(I; x = xr, y = yr, z = zr, level = level, cap = true, padValue=1e8)      

        fig = Figure(size=(1600,800))
        ax1 = AxisGeom(fig[1, 1], title = "gradient vectors")
        hp1 = meshplot!(ax1, F, V, color=(:white, 0.9), transparency=true, strokewidth=0.0)        
            # hpa = dirplot(ax1,V,N; color=:blue,linewidth=3,scaleval=1.0,style=styleSet[i])
        arrows3d!(ax1, P0, P, lengthscale = 0.5, color = C)

        screen = display(GLMakie.Screen(), fig)
        GLMakie.set_title!(screen, "testCase = $testCase")
    end
end