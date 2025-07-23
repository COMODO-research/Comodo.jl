using Comodo
using Comodo.GLMakie

#=
This demo shows the use of `pointspacingmax` function to obtain the largest 
point spacing of the input curve, face, or faces. 
=#

for testCase = 1:4
    if testCase == 1 # Curve not closed
        V = Point{3,Float64}[[0,0,0],[2.0,0,0],[1,1,0],[0,1,0]]        
        d_true = 2.0
        d = pointspacingmax(V; close_loop=false)
    elseif testCase == 2 # Curve closed        
        w = 2.5
        F, V = square(w)  
        V[1] = Point{3,Float64}(w/2.0,w,0.0)       
        d_true = w+w/2.0
        d = pointspacingmax(V; close_loop=true)
    elseif testCase == 3 # Single face
        w = 3.0
        f, V = square(w)  
        V[1] = Point{3,Float64}(w/2.0,w,0.0)       
        d_true = w+w/2.0
        d = pointspacingmax(f,V)
    elseif testCase == 4 # Multiple faces
        w = 3.0 # width
        F, V = cube(w*sqrt(3)/2.0) 
        for i in F[1]
            V[i] = V[i].*2.0
        end
        d_true = w*2.0
        close_loop = true
        d = pointspacingmax(F,V)
    end    
    println("Max length: " * string(d) * ", theoretical/true: " * string(d_true))
end
