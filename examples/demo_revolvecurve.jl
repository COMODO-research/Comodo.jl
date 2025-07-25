using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

GLMakie.closeall()

for testCase = 1:5
    if testCase == 1
        θ = 1.0*pi
        r = 1.0
        nc = 15
        t = range(0,2*π,nc)
        Vc = [Point3{Float64}(2.0+tt,0.0,sin(tt)) for tt ∈ t]
        Vc = evenly_sample(Vc, nc)
        n = Vec{3, Float64}(0.0,0.0,1.0)
        num_steps = 31
        periodicity = (false,false)
    elseif testCase == 2
        θ = 1.5*pi
        r = 1.0
        nc = 15
        t = range(0,2*π,nc)
        Vc = [Point3{Float64}(2.0+tt,0.0,sin(tt)) for tt ∈ t]
        Vc = evenly_sample(Vc, nc)
        n = normalizevector(Vec{3, Float64}(0.0,-1.0,1.0))
        num_steps = 31
        periodicity = (false,false)
    elseif testCase == 3
        θ = 1.0*pi
        r = 2.0
        nc = 16
        t = range(2.0*π-(2.0*π/nc),0,nc)
        Vc = [Point3{Float64}(4.0+r*cos(tt),0.0,r*sin(tt)) for tt ∈ t]    
        n = Vec{3, Float64}(0.0,0.0,1.0)
        num_steps = 25
        periodicity = (false,false)
    elseif testCase == 4
        θ = 1.0*pi
        r = 2.0
        nc = 16
        t = range(2.0*π-(2.0*π/nc),0,nc)
        Vc = [Point3{Float64}(4.0+r*cos(tt),0.0,r*sin(tt)) for tt ∈ t]     
        n = Vec{3, Float64}(0.0,0.0,1.0)
        num_steps = 25
        periodicity = (true,false) 
    elseif testCase == 5        
        r = 2.0
        nc = 16
        num_steps = 50
        θ = 2.0*pi-(2.0*pi/num_steps)
        t = range(2.0*π-(2.0*π/nc),0,nc)
        Vc = [Point3{Float64}(4.0+r*cos(tt),0.0,r*sin(tt)) for tt ∈ t]    
        n = Vec{3, Float64}(0.0,0.0,1.0)    
        periodicity = (true,true)    
    end

    face_type1 =:quad
    face_type2 =:forwardslash
    face_type3 =:backslash
    face_type4 =:tri
    face_type5 =:tri_even

    direction1 = :positive
    direction2 = :both
    direction3 = :negative
    direction4 = :positive
    direction5 = :positive

    F1,V1 = revolvecurve(Vc; extent=θ, direction=direction1, n=n ,num_steps=num_steps, periodicity=periodicity,face_type=face_type1)
    F2,V2 = revolvecurve(Vc; extent=θ, direction=direction2, n=n, num_steps=num_steps, periodicity=periodicity,face_type=face_type2)
    F3,V3 = revolvecurve(Vc; extent=θ, direction=direction3, n=n, num_steps=num_steps, periodicity=periodicity,face_type=face_type3)
    F4,V4 = revolvecurve(Vc; extent=θ, direction=direction4, n=n, num_steps=num_steps, periodicity=periodicity,face_type=face_type4)
    F5,V5 = revolvecurve(Vc; extent=θ, direction=direction5, n=n, num_steps=num_steps, periodicity=periodicity,face_type=face_type5)
    F6,V6 = revolvecurve(Vc; extent=2*pi-(2*pi)/num_steps, direction=direction5, n=n, num_steps=num_steps, periodicity=(false,true),face_type=face_type5)

    ## Visualization
    markersize = 8 
    linewidth = 3

    fig = Figure(size=(1800,800))

    ax1 = AxisGeom(fig[1, 1], title = """Revolved $θ rad, direction=$direction1, periodicity=$periodicity, face_type=$face_type1""")
    hp1 = lines!(ax1,Vc,color=:red,linewidth=linewidth, transparency=true, depth_shift=-1.0f-3)
    hp2 = scatter!(ax1,Vc,markersize=markersize,color=:red)
    hp3 = meshplot!(ax1, F1, V1)
    arrows!(ax1,[Point{3,Float64}(0.0,0.0,0.0)], [7*n], fxaa=true, color=:green,arrowsize = 1)

    ax2 = AxisGeom(fig[1, 2], title = """Revolved $θ rad, direction=$direction2, periodicity=$periodicity, face_type=$face_type2""")
    hp1 = lines!(ax2,Vc,color=:red,linewidth=linewidth, transparency=true, depth_shift=-1.0f-3)
    hp2 = scatter!(ax2,Vc,markersize=markersize,color=:red)
    hp3 = meshplot!(ax2, F2, V2)
    arrows!(ax2,[Point{3,Float64}(0.0,0.0,0.0)], [7*n], fxaa=true, color=:green,arrowsize = 1)

    ax3 = AxisGeom(fig[1, 3], title = """Revolved $θ rad, direction=$direction3, periodicity=$periodicity, face_type=$face_type3""")
    hp1 = lines!(ax3,Vc,color=:red,linewidth=linewidth, transparency=true, depth_shift=-1.0f-3)
    hp2 = scatter!(ax3,Vc,markersize=markersize,color=:red)
    hp3 = meshplot!(ax3, F3, V3)
    arrows!(ax3,[Point{3,Float64}(0.0,0.0,0.0)], [7*n], fxaa=true, color=:green,arrowsize = 1)

    ax4 = AxisGeom(fig[2, 1], title = """Revolved $θ rad, direction=$direction4, periodicity=$periodicity, face_type=$face_type4""")
    hp1 = lines!(ax4,Vc,color=:red,linewidth=linewidth, transparency=true, depth_shift=-1.0f-3)
    hp2 = scatter!(ax4,Vc,markersize=markersize,color=:red)
    hp3 = meshplot!(ax4, F4, V4)
    arrows!(ax4,[Point{3,Float64}(0.0,0.0,0.0)], [7*n], fxaa=true, color=:green,arrowsize = 1)

    ax5 = AxisGeom(fig[2, 2], title = """Revolved $θ rad, direction=$direction5, periodicity=$periodicity, face_type=$face_type5""")
    hp1 = lines!(ax5,Vc,color=:red,linewidth=linewidth, transparency=true, depth_shift=-1.0f-3)
    hp2 = scatter!(ax5,Vc,markersize=markersize,color=:red)
    hp3 = meshplot!(ax5, F5, V5)
    arrows!(ax5,[Point{3,Float64}(0.0,0.0,0.0)], [7*n], fxaa=true, color=:green,arrowsize = 1)

    ax6 = AxisGeom(fig[2, 3], title = """Revolved $θ rad, direction=$direction5, periodicity=$periodicity, face_type=$face_type5""")
    hp1 = lines!(ax6,Vc,color=:red,linewidth=linewidth, transparency=true, depth_shift=-1.0f-3)
    hp2 = scatter!(ax6,Vc,markersize=markersize,color=:red)
    hp3 = meshplot!(ax6, F6, V6)
    arrows!(ax6,[Point{3,Float64}(0.0,0.0,0.0)], [7*n], fxaa=true, color=:green,arrowsize = 1)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end