using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics

#=
This demo shows the use of hermiteSegment to set up a Hermite spline using a set
of two control points and two velocity vectors. "Under the hood" this function
simply uses a Bezier spline with the 4 control points: 
P = [p1, p1+v1/norm(v1).0, p2-v2/norm(v2).0, p2]. 
=#

n = 40
p1 = Point{3,Float64}(0.0, 0.0, 0.0)
p2 = Point{3,Float64}(1.0, 1.0, 0.0)
v1 = Point{3,Float64}(0.0, 1.0, 0.0)
v2 = Point{3,Float64}(1.0, .0, 0.0)
V = hermiteSegment(n, p1, v1, p2, v2)

# Visualization
vMax = 6.0
fig = Figure(size = (800,800))
ax = AxisGeom(fig[1, 1], azimuth=-pi/2, elevation=pi/2, limits=(-0.5,1.5,-0.5,1.5,0,1))

hp1 = scatter!(ax, [p1,p2], markersize=25,color=:black)
hp2 = lines!(ax, V,linewidth=5,color=:red)
hp3 = lines!(ax,[p1,p1+v1/norm(v1)/2.0]; color=norm(v1), linewidth=3, colormap=:viridis, colorrange=(0.0,vMax))
hp4 = lines!(ax,[p2,p2+v2/norm(v2)/2.0]; color=norm(v2), linewidth=3, colormap=:viridis, colorrange=(0.0,vMax))
# arrows3d!(ax,p1, v1, shaftcolor = :green, tipcolor = :green, align = :center)
Legend(fig[1, 2],[hp1, hp2, hp3],["Control points","Hermite (Bezier) spline", "Velocity vector"])

function updatePlot(V,p1,p2,v1,v2)
    hp2[1] = V   
    hp3[1] = [p1,p1+v1/norm(v1)/2.0]    
    hp4[1] = [p2,p2+v2/norm(v2)/2.0]
    hp3.color = norm(v1)
    hp4.color = norm(v2)
end

stepRange = range(0.0,vMax,100)
hSlider1 = Slider(fig[2, :], range = stepRange, startvalue = 1.0, linewidth=30)
on(hSlider1.value) do v     
    a1 = hSlider3.value[]    
    v1m = v
    v1 = Point{3,Float64}(v1m*cosd(a1), v1m*sind(a1), 0.0)

    a2 = hSlider4.value[]    
    v2m = hSlider2.value[]
    v2 = Point{3,Float64}(v2m*cosd(a2), v2m*sind(a2), 0.0)

    V = hermiteSegment(n, p1, v1, p2, v2)
    updatePlot(V,p1,p2,v1,v2)    
end

hSlider2 = Slider(fig[3, :], range = stepRange, startvalue = 1.0, linewidth=30)
on(hSlider2.value) do v     
    a1 = hSlider3.value[]    
    v1m = hSlider1.value[]
    v1 = Point{3,Float64}(v1m*cosd(a1), v1m*sind(a1), 0.0)

    a2 = hSlider4.value[]    
    v2m = v
    v2 = Point{3,Float64}(v2m*cosd(a2), v2m*sind(a2), 0.0)

    V = hermiteSegment(n, p1, v1, p2, v2)
    updatePlot(V,p1,p2,v1,v2)
end

stepRange = range(-360,360.0,200)
hSlider3 = Slider(fig[4, :], range = stepRange, startvalue = 90, linewidth=30)
on(hSlider3.value) do a     
    a1 = a    
    v1m = hSlider1.value[]
    v1 = Point{3,Float64}(v1m*cosd(a1), v1m*sind(a1), 0.0)

    a2 = hSlider4.value[]    
    v2m = hSlider2.value[]
    v2 = Point{3,Float64}(v2m*cosd(a2), v2m*sind(a2), 0.0)

    V = hermiteSegment(n, p1, v1, p2, v2)
    updatePlot(V,p1,p2,v1,v2)
end

hSlider4 = Slider(fig[5, :], range = stepRange, startvalue = 0, linewidth=30)
on(hSlider4.value) do a     
    a1 = hSlider3.value[]    
    v1m = hSlider1.value[]
    v1 = Point{3,Float64}(v1m*cosd(a1), v1m*sind(a1), 0.0)

    a2 = a
    v2m = hSlider2.value[]
    v2 = Point{3,Float64}(v2m*cosd(a2), v2m*sind(a2), 0.0)

    V = hermiteSegment(n, p1, v1, p2, v2)
    updatePlot(V,p1,p2,v1,v2)
end

fig
