using Comodo
using Comodo.GLMakie
using Comodo.Rotations
#=
This demo shows the use of `fitplane` to fit a plane of the for z=Ax+By+C to a set
of points `V`. 
=#

GLMakie.closeall() # Close all figures 

## Create example plane
plateDim = [4.0, 8.0]
plateElem = [12, 12]
sd = 0.25 # Standard deviation for noise in the below demos
for testCase = 1:3
    if testCase==1 # x=Ay+Bz+C
        F, V, _, _ = quadplate(plateDim, plateElem; orientation=:up)

        # Define a rotation tensor using Euler angles
        Q = RotXYZ(0.0, 0.5*π, 0.0)

        # Rotate the coordinates
        V = [Point{3, Float64}(Q*v) for v ∈ V] 

        A_true = 0.25
        B_true = -0.3
        C_true = -5.0

        Vt = [Point{3, Float64}((A_true*p[2]+B_true*p[3]+C_true), p[2], p[3]) for p in V]
        Vtn = [Point{3, Float64}(p[1], p[2], p[3] + sd*randn()) for p in Vt]

        ## Now fit plane of form x=Ay+Bz+C
        A, B, C = fitplane(Vtn; direction=1)

        Vf = [Point{3, Float64}(A*p[2]+B*p[3]+C, p[2], p[3]) for p in V]        
        f = "x=Ay+Bz+C"
    elseif testCase==2 # y=Ax+Bz+C
        F, V, _, _ = quadplate(plateDim, plateElem; orientation=:up)

        # Define a rotation tensor using Euler angles
        Q = RotXYZ(0.5*π, 0.0, 0.0)

        # Rotate the coordinates
        V = [Point{3, Float64}(Q*v) for v ∈ V] 

        A_true = 0.25
        B_true = -0.3
        C_true = -5.0
        
        Vt = [Point{3, Float64}(p[1], (A_true*p[1]+B_true*p[3]+C_true), p[3]) for p in V]
        Vtn = [Point{3, Float64}(p[1], p[2], p[3] + sd*randn()) for p in Vt]

        ## Now fit plane of form y=Ax+Bz+C
        A, B, C = fitplane(Vtn; direction=2)

        Vf = [Point{3, Float64}(p[1], A*p[1]+B*p[3]+C, p[3]) for p in V]  
        f = "y=Ax+Bz+C"
    elseif testCase==3 # z=Ax+By+C        
        F, V, _, _ = quadplate(plateDim, plateElem; orientation=:up)

        A_true = 0.25
        B_true = -0.3
        C_true = -5.0

        Vt = [Point{3, Float64}(p[1], p[2], (A_true*p[1]+B_true*p[2]+C_true)) for p in V]
        Vtn = [Point{3, Float64}(p[1], p[2], p[3] + sd*randn()) for p in Vt]

        ## Now fit plane of form z=Ax+By+C
        A, B, C = fitplane(Vtn; direction=3)

        Vf = [Point{3, Float64}(p[1], p[2], A*p[1]+B*p[2]+C) for p in V]
        f = "z=Ax+By+C"
    end

    println("True: ", [A_true, B_true, C_true])
    println("Fit: ", [A, B, C])

    ## Visualisation
    fig = Figure(size=(1600,800))

    ax1 = AxisGeom(fig[1, 1], title = "Fitted plane of form: $f, A=$A, B=$B, C=$C")        
    scatter!(ax1, Vtn, markersize=15, color=:black)
    hp2 = meshplot!(ax1, F, Vt; strokewidth=0.0, color=(:green, 0.5), transparency=true)

    hp2 = meshplot!(ax1, F, Vf; strokewidth=0.0, color=(:red, 0.5), transparency=true)

    screen = display(GLMakie.Screen(), fig)
end