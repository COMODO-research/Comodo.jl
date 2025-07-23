using Comodo
using Comodo.GLMakie

#=
This demo shows the use of `pointspacingmean` function to obtain the mean point 
spacing of the input curve, face, or faces. 
=#

for testCase = 1:4
    if testCase == 1 # Curve not closed
        V = Point{3,Float64}[[0,0,0],[1.0,0,0],[1,1,0],[0,1,0]]        
        d_true = 1
        d = pointspacingmean(V; close_loop=false)
    elseif testCase == 2 # Curve closed
        d_true = 2.0 # width
        F, V = square(d_true)         
        d = pointspacingmean(V; close_loop=true)
    elseif testCase == 3 # Single face
        d_true = 2.0 # width
        f, V = square(d_true)         
        d = pointspacingmean(f,V)
    elseif testCase == 4 # Multiple faces
        d_true = 3.0 # width
        F, V = cube(d_true*sqrt(3)/2.0) 
        close_loop = true
        d = pointspacingmean(F,V)
    end    
    println("Mean length: " * string(d) * ", theoretical/true: " * string(d_true))
end