using Comodo
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.GLMakie
using Comodo.Distributions
using Comodo.Rotations 
using Random

GLMakie.closeall()

function vec_cone(m::Point{3,Float64}, r=1.0, β=π/4, N=25)
    nz = Point{3,Float64}(0.0, 0.0, 1.0)
    if m != nz
        Qm = rotation_between(nz,m)         
        rotate_bool = true
    else
        rotate_bool = false
    end
    V = Vector{Point{3,Float64}}(undef, N+1)    
    Q = RotX(β)
    n = Q*nz    
    a = -(2*pi)/(N)
    if rotate_bool 
        V[1] = Qm*r*n
    else
        V[1] = r*n
    end
    V[end] = Point{3,Float64}(0.0, 0.0, 0.0)
    for i in 2:N        
        Q = RotZ((i-1)*a)        
        if rotate_bool 
            Q=Qm*Q
        end
        V[i] = r*Q*n        
    end    
    F = [TriangleFace(i, mod1(i+1,N), N+1) for i in 1:N] # Construct faces
    return F, V
end

r = 1.0
N = 1500
β = π/6

θn = 0.25*π
ϕn = -0.3*π
n = Point{3,Float64}(cos(ϕn)*sin(θn), sin(ϕn)*sin(θn), cos(θn))
P = rand_onsphere_cone(r, β, N, n)

# Visualisation

fig = Figure(size=(800,800))
ax1 = AxisGeom(fig[1, 1], title="random points on sphere in cone")

F, V = geosphere(3, r)
Fc1, Vc1 = vec_cone( n, r, β, 50)
hp3 = meshplot!(ax1, F, V, strokewidth=0.0, color=(:white,0.5), transparency=true)
meshplot!(ax1, Fc1, Vc1, strokewidth=0.0, color=(:red,0.5), transparency=true)
scatter!(ax1, P, markersize=5, color=:black)
dirplot(ax1, Point{3,Float64}(0.0, 0.0, 0.0), n, linewidth=4, color=:black)
screen = display(GLMakie.Screen(), fig)
