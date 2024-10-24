using Comodo
using GLMakie
using GeometryBasics
using LinearAlgebra

# Example data 


convert_method = :angle

testCase = 1
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

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original")
poly!(ax1,GeometryBasics.Mesh(Vs,Fs), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Forward slash converted")
poly!(ax2,GeometryBasics.Mesh(Vs1,Fs1), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax3 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Back slash converted")
poly!(ax3,GeometryBasics.Mesh(Vs2,Fs2), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax4 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Angle based conversion")
poly!(ax4,GeometryBasics.Mesh(Vs3,Fs3), strokewidth=3,color=:white,shading=FastShading,transparency=false)

fig