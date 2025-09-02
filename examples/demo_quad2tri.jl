using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

GLMakie.closeall()

for testCase = 1:3
    if testCase ==1
        r = 1.0 # Radius
        n = 3 # Number of refinement steps from cube
        F,V = quadsphere(n,r)
    elseif testCase ==2
        F,V = quadplate([5,5],[5,5])
        for i ∈ eachindex(V)
            if V[i][2]>1
                V[i]=[V[i][1] + V[i][2]*0.5 V[i][2] V[i][3]]
            elseif V[i][2]<-1
                V[i]=[V[i][1] + V[i][2]*-0.5 V[i][2] V[i][3]]
            end
        end
    elseif testCase ==3
        F,V = cube(1.0)
        
        # Build deformation gradient tensor to induce shear with known angles
        fd = zeros(3,3)
        for i=1:3
            fd[i,i]=1.0
        end
        a = pi/6 # "45 degree shear"  
        fd[1,2] = tan(a) 
        V = [eltype(V)(fd*v) for v ∈ V] # Shear the cube

        # Build deformation gradient tensor to induce shear with known angles
        fd = zeros(3,3)
        for i=1:3
            fd[i,i]=1.0
        end
        a = pi/6 # "45 degree shear"  
        fd[3,1] = tan(a) 
        
        V = [eltype(V)(fd*v) for v ∈ V] # Shear the cube

    end

    Ft1 = quad2tri(F,V; convert_method = :forward)
    Ft2 = quad2tri(F,V; convert_method = :backward)
    Ft3 = quad2tri(F,V; convert_method = :angle)

    ## Visualization

    Fs,Vs = separate_vertices(F,V)
    Fs1,Vs1 = separate_vertices(Ft1,V)
    Fs2,Vs2 = separate_vertices(Ft2,V)
    Fs3,Vs3 = separate_vertices(Ft3,V)

    fig = Figure(size=(800,800))

    ax1 = AxisGeom(fig[1, 1], title = "Original")
    meshplot!(ax1, Fs, Vs, strokewidth=3)

    ax2 = AxisGeom(fig[1, 2], title = "Forward slash converted")
    meshplot!(ax2,Fs1, Vs1, strokewidth=3)

    ax3 = AxisGeom(fig[2, 1], title = "Back slash converted")
    meshplot!(ax3,Fs2, Vs2, strokewidth=3)

    ax4 = AxisGeom(fig[2, 2], title = "Angle based conversion")
    meshplot!(ax4,Fs3, Vs3, strokewidth=3)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end