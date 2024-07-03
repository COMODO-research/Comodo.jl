using Comodo
using GLMakie
using GeometryBasics
using LinearAlgebra
using Rotations
using Random


# Set seed so demo performs the same each time
Random.seed!(1)

testCase = 11

if testCase == 1
    V = Point{3,Float64}[ [0.0,0.0,0.0], [10.0,0.0,0.0]]
    close_loop = false
elseif testCase == 2
    d = 10.0
    V = Point{3,Float64}[ [d,0.0,0.0], [d,d,0.0], [0.0,d,0.0]]
    close_loop = false
elseif testCase == 3
    d = 10.0 
    V = Point{3,Float64}[ [-d,-d,0.0], [d,-d,0.0], [d,d,0.0], [-d,d,0.0]]
    close_loop = false
elseif testCase == 4
    d = 10.0 
    V = Point{3,Float64}[ [-d,-d,0.0], [d,-d,0.0], [d,d,0.0], [-d,d,0.0]]
    close_loop = true
elseif testCase == 5
    V = Point{3,Float64}[ [0.0,0.0,0.0], [1.0,0.0,0.0], [1.0,2.0,0.0], [1.0,2.0,1.0]]
    close_loop = false
elseif testCase == 6 
    d = 10.0 
    nc = 12
    V = circlepoints(d,nc)  
    close_loop = false
elseif testCase == 7 
    V = Point{3,Float64}[ [-2.0,0.0,0.0], [-1.0,0.0,0.0], [0.0,0.0,0.0], [1.0,-1.0,1.0], [2.0,1.0,0.0], [3.0,2.0,8.0], [5.0,0.0,0.0], [7.0,0.0,0.0]]
    close_loop = false
elseif testCase == 8 
    close_loop = false
    V = 6*[Point{3,Float64}(i,6*rand(1)[1],6*rand(1)[1]) for i in 1:1:20]
elseif testCase == 9
    V = circlepoints(20,24; dir=:acw)
    V = [Point{3,Float64}(v[1],v[2],v[3]+20*rand(1)[1]) for v in V]
    V = [v+12*rand(3) for v in V]
    close_loop = true
elseif testCase == 10
    V = circlepoints(20,24; dir=:acw)
    for i = 1:2:length(V)
        v = V[i]
        V[i] = Point{3,Float64}(v[1],v[2],v[3]+20.0)
    end
    close_loop = true
elseif testCase == 11
    d = 10.0 
    V = Point{3,Float64}[ [-d,-d,0.0], [d,-d,0.0], [d,d,0.0], [-d,d,0.0]]
    close_loop = true
end


rMax = [2.0,3.0,0.0,5.0] #nothing

VC = filletcurve(V; rMax=rMax,  constrain_method = :max, n=25, close_loop = close_loop, eps_level = 1e-6)
# VC2,_ = evenly_sample(VC,100)

# Visualisation
fig = Figure(size=(1000,1000))
ax1 = Axis3(fig[1, 1],aspect = :data,title="Input curve")

hp11 = lines!(ax1, V,linewidth=2,color=:black)
if close_loop == true
    hp12 = lines!(ax1, [V[1],V[length(V)]],linewidth=2,color=:black)
end
hp2 = scatter!(ax1, V,markersize=15,color=:black)

# scatter!(ax1, V[1],markersize=25,color=:yellow)
# scatter!(ax1, V[end],markersize=25,color=:red)

indPlot = collect(1:length(VC))
if close_loop == true
    push!(indPlot,1)
end
hp2 = lines!(ax1, VC[indPlot],linewidth=6,color=:blue)

stepRange1 = range(0,10,500)
hSlider1 = Slider(fig[2, 1], range = stepRange1, startvalue = 0,linewidth=30)

on(hSlider1.value) do stepIndex1
    VC = filletcurve(V; rMax=stepIndex1,  constrain_method = :max, n=25, close_loop = close_loop, eps_level = 1e-6)
    indPlot = collect(1:length(VC))
    if close_loop == true
        push!(indPlot,1)
    end
    hp2[1] = VC[indPlot]
end

fig