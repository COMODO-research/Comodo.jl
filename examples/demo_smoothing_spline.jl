using Comodo
using Comodo.GeometryBasics
using Comodo.GLMakie
using Comodo.BSplineKit
using Printf

#=
This demo uses BSplineKit to fit a cubic smoothing spline to a parameterised 
curve. In this case the data represents noisy radius data for a circular curve. 
Two cases are implemented here, for testCase = 1 the 3D coordinates are 
parameterised with respect to curve length. If instead textCase = 2 then the 
curve is parameterised with respect to angle and the radii are fitted. In the 
latter case the smoothest result is therefore a circle fit. 
=#

n = 80
t = collect(range(0.0,(2*pi)-(2*pi)/n,n))
r = 1.0 .+ 0.1.*randn(n) # Noisy radius data
V = Vector{Point{3,Float64}}(undef,n)
for (i,tt) in enumerate(t)
    V[i] = Point{3,Float64}(r[i]*cos(tt),r[i]*sin(tt),0.0)
end

λ = 0.0 # smoothing parameter

m = 1000

GLMakie.closeall()

for testCase = 1:2
    if testCase == 1
        L = curve_length(V; close_loop=true)
        S = BSplineKit.fit(L[1:end-1], V, λ, BSplineKit.Periodic(maximum(L)))
        L_fit = range(0.0,maximum(L),m)
        V_fit = S.(L_fit)
    elseif testCase == 2        
        S = BSplineKit.fit(t, r, λ, BSplineKit.Periodic(2.0*pi))
        t_fit = range(0.0,2*pi,m)
        r_fit = S.(t_fit)
        V_fit = Vector{Point{3,Float64}}(undef,m)
        for (i,tt) in enumerate(t_fit)
            V_fit[i] = Point{3,Float64}(r_fit[i]*cos(tt),r_fit[i]*sin(tt),0.0)
        end
    end

    # Visualization
    fig = Figure(size = (800,800))    
    ax = AxisGeom(fig[1, 1])

    hp1 = scatter!(ax, V,markersize=15,color=:red)
    hp2 = lines!(ax, V_fit,linewidth=4,color=:black)

    Legend(fig[1, :2],[hp1,hp2],["Raw","Smoothing spline"])

    stepRange = [0.0,0.0,0.000025,0.00005,0.0001,0.00015,0.0002,0.0005,0.0015,0.003,0.02,0.05,0.1,0.5,1];#range(0,0.5,1000)
    hSlider = Slider(fig[2,:], range = stepRange, startvalue = 0, linewidth=30)

    on(hSlider.value) do λ
        if testCase == 1
            S = BSplineKit.fit(L[1:end-1], V, λ, BSplineKit.Periodic(maximum(L)))    
            V_fit = S.(L_fit)
        elseif testCase == 2  
            S = BSplineKit.fit(t, r, λ, BSplineKit.Periodic(2*pi))
            r_fit = S.(t_fit)
            V_fit = Vector{Point{3,Float64}}(undef,m)
            for (i,tt) in enumerate(t_fit)
                V_fit[i] = Point{3,Float64}(r_fit[i]*cos(tt),r_fit[i]*sin(tt),0.0)
            end
        end

        hp2[1] = V_fit    
        
        ax.title =  @sprintf "Smoothing spline λ = %f" λ  
    end
    slidercontrol(hSlider,ax)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end