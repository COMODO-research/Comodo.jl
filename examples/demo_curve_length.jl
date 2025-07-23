using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.BSplineKit
using Comodo.QuadGK: quadgk # For numerical integration

for testCase = 1:3
    # Define input curve
    if testCase == 1
        V = Point{3,Float64}[[0,0,0],[1.0,0,0],[1,1,0],[0,1,0]]
        close_loop = false
        d_true = 3
    elseif testCase == 2
        V = Point{3,Float64}[[0,0,0],[1.0,0,0],[1,1,0],[0,1,0]]
        close_loop = true
        d_true = 4
    elseif testCase == 3
        n = 250
        r = 3.5
        V = circlepoints(r,n; dir=:acw)
        close_loop = true
        d_true = 2*pi*r
    end
    L = curve_length(V; close_loop = close_loop) # Allong curve distance   
    d = curve_length(V; close_loop = close_loop, total=true) # Total distance  
    println("Max. length: " * string(d) * ", theoretical/true: " * string(d_true))
end