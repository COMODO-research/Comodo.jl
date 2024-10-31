using Comodo
using GLMakie
using GeometryBasics

#=
This demo shows the use of `hexbox` to generate a hexahedral mesh for a 3D box
domain. 
=#

testCase = 5

if testCase == 1 # Batman curve
    n = 120
    V = batman(n; symmetric=true)
    pointSpacing = pointspacingmean(V)

    VT = (V,)
    R = ([1],)
    P = (pointSpacing)
elseif testCase == 2 
    n = 50
    r = 2.0
    V = circlepoints(r,n;dir=:acw) 
    pointSpacing = pointspacingmean(V)

    VT = (V,)
    R = ([1],)
    P = (pointSpacing)
elseif testCase == 3 
    n = 100
    r = 1.0
    rFun(t) = r + 0.5.*sin(6*t)
    V = circlepoints(rFun,n; dir=:acw)    
    V = evenly_sample(V, n)
    pointSpacing = pointspacingmean(V)

    VT = (V,)
    R = ([1],)
    P = (pointSpacing)
elseif testCase == 4
    rFun1(t) = 12.0 + 3.0.*sin(5*t)
    n1 = 120
    V1 = circlepoints(rFun1,n1)
    V1 = evenly_sample(V1,n1)

    rFun2(t) = 10.0 + 2.0.*sin(5*t)
    n2 = 100
    V2 = circlepoints(rFun2,n2)
    V2 = evenly_sample(V2,n2)

    rFun3(t) = 2.0 + 0.5 *sin(5*t)
    n3 = 50
    Vp = circlepoints(rFun3,n3)
    Vp = evenly_sample(Vp,n3)

    V3 = [Point{3,Float64}(v[1],v[2]+6,v[3]) for v in Vp]
    V4 = [Point{3,Float64}(v[1]-5,v[2]+1,v[3]) for v in Vp]
    V5 = [Point{3,Float64}(v[1]+5,v[2]+1,v[3]) for v in Vp]
    V6 = [Point{3,Float64}(v[1],v[2]-4,v[3]) for v in Vp]
    
    VT = (V1,V2,V3,V4,V5,V6) # Curves
    R = ([1,2],[2,3,4,5,6],[4],[5]) # Regions 
    P = (1,0.6,0.2,0.3)  # Point spacings

elseif testCase == 5    
    n1 = 120
    r1 = 20.0
    V1 = circlepoints(r1,n1)
    
    n2 = 100
    r2 = 12.0
    V2 = circlepoints(r2,n2)
    V2 = [Point{3,Float64}(v[1]+6,v[2],v[3]) for v in V2]
        
    n3 = 30
    r3 = 4
    V3 = circlepoints(r3,n3)
    V3 = [Point{3,Float64}(v[1]-14,v[2],v[3]) for v in V3]
        
    n4 = 50
    r4 = 7
    V4 = circlepoints(r4,n4)
    V4 = [Point{3,Float64}(v[1]+9,v[2],v[3]) for v in V4]

    n5 = 40
    r5 = 3
    V5 = circlepoints(r5,n5)
    V5 = [Point{3,Float64}(v[1]-2,v[2],v[3]) for v in V5]

    VT = (V1,V2,V3,V4,V5) # Curves
    R = ([1,2,3],[2,4,5],[5]) # Regions 
    P = (1,0.75,0.5)  # Point spacings
end

F,V,C = regiontrimesh(VT,R,P)


# Visualisation

Fp,Vp = separate_vertices(F,V) # Give each face its own point set 
Cp = simplex2vertexdata(Fp,C) # Convert face color data to vertex color data 


fig = Figure(size=(1000,1000))
ax1 = Axis3(fig[1, 1],aspect = :data,title="Multi-region meshing",azimuth=-pi/2,elevation=pi/2)

hp4 = poly!(ax1,GeometryBasics.Mesh(Vp,Fp), strokewidth=1,color=Cp, strokecolor=:black, shading = FastShading, transparency=false,colormap=Makie.Categorical(Makie.Reverse(:Spectral)))

fig
