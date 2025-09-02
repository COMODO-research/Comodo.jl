using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Printf

function truncatedntrapezohedron(n,r=1.0,f=nothing)
    m = 2*n # Number of faces
    R = 0.5*csc(pi/n)*r # Radius scaled to obtain desired midradius
    zc = r * sqrt(4-sec(pi/m)^2)/(4+8*cos(pi/n)) # z-offset of equatorial points
    zt = r * 1/4*cos(pi/m)*cot(pi/m)*csc(3*pi/m)*sqrt(4-sec(pi/m)^2) # z-offset of poles

    # Create equatorial points
    T = circlerange(m; dir=:acw) # Range of angles for points in circle
    V = Vector{Point{3,Float64}}(undef,2*m)
    for (i,t) in enumerate(T)
        if iseven(i)
            s = 1.0
        else
            s = -1.0
        end        
        V[i] = Point{3, Float64}(R*cos(t),R*sin(t),s*zc)     
    end
    d1 = norm(V[1]-V[2])
    d2 = norm(V[1]-Point{3, Float64}(0.0,0.0,-zt))

    if isnothing(f)
        f=1.0-(d1/d2)
    end

    for i=1:n
        j = 2 + (i-1)*2        
        V[m+i  ] = f*V[j-1] + (1.0-f) * Point{3, Float64}(0.0,0.0,-zt)
        V[m+i+n] = f*V[j] + (1.0-f) * Point{3, Float64}(0.0,0.0,zt)     
    end

    # Construct faces
    F1 =  Vector{NgonFace{5,Int}}(undef,m)
    for i = 1:n # Use n steps to create m=2*n faces
        j = 1 + 2*(i-1) # First point index for bottom pentahedra          
        F1[i  ] = NgonFace{5,Int}(mod1(j+2,m), mod1(j+1,m), j, m+i, m+mod1(i+1,n))  # Bottom pentahedra 
        F1[i+n] = NgonFace{5,Int}(mod1(j+1,m), mod1(j+2,m), mod1(j+3,m),  m+n+mod1(i+1,n), m+n+i)  # Top pentahedra     
    end

    F2 =  Vector{NgonFace{n,Int}}(undef,2)
    F2[1] = NgonFace{n,Int}(m+n:-1:m+1)  # Bottom n-hedra 
    F2[2] = NgonFace{n,Int}(m+n+1:2*m)  # Bottom n-hedra 

    F = (F1,F2)

    return F,V
end

n = 12
r = 1.0
f = nothing

F,V = truncatedntrapezohedron(n,r,f)

## Visualize mesh
GLMakie.closeall()

markersize = 25
strokewidth = 2 
strokecolor = :black

fig = Figure(size = (1200,1200))
ax1 = AxisGeom(fig[1, 1], title = "truncated n-trapezohedron")

hp1 = meshplot!(ax1, GeometryBasics.Mesh(V,F[1]))
hp3 = scatter!(ax1, V,markersize=markersize,color=:red)

stepRange = 3:12
hSlider = Slider(fig[2,1], range = stepRange, startvalue = 0,linewidth=30)

on(hSlider.value) do n 
    F,V = truncatedntrapezohedron(n,r,f)
    
    M1 = GeometryBasics.Mesh(V,F[1])
    # M2 = GeometryBasics.Mesh(V,F[2])

    hp1[1] = M1    
    # hp2[1] = M2
    hp3[1] = V    
    ax1.title =  @sprintf "truncated n-trapezohedron n=%i" n  
end
slidercontrol(hSlider,ax1)

fig