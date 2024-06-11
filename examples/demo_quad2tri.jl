using Comodo
using GLMakie
using GeometryBasics
using LinearAlgebra

# Example data 


convert_method = :angle

testCase = 3
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
    M = cube(1.0)
    F = faces(M)
    V = coordinates(M)

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

Ft = quad2tri(F,V; convert_method = convert_method)


## Visualization
fig = Figure(size=(800,800))

ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Original")
poly!(ax1,GeometryBasics.Mesh(V,F), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax2 = Axis3(fig[1, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Forward slash converted")
poly!(ax2,GeometryBasics.Mesh(V,quad2tri(F,V; convert_method = :forward)), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax3 = Axis3(fig[2, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Back slash converted")
poly!(ax3,GeometryBasics.Mesh(V,quad2tri(F,V; convert_method = :backward)), strokewidth=3,color=:white,shading=FastShading,transparency=false)

ax4 = Axis3(fig[2, 2], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Angle based conversion")
poly!(ax4,GeometryBasics.Mesh(V,quad2tri(F,V; convert_method = :angle)), strokewidth=3,color=:white,shading=FastShading,transparency=false)

fig