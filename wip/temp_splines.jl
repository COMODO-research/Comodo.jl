using GLMakie
# using BenchmarkTools
# using Makie.GeometryBasics
# using ColorSchemes

n = 25

p = [0 0 0; 1 1 1; 2 -1 1; 3 0 0;] # Control points

function BezierFactorial(p,n)
    t   = collect(range(0,1.0,n))
    N   = size(p,1)
    nn  = collect(0:1:N-1)'
    f   = factorial.(nn)
    s   = factorial(N-1)./( f.*reverse(f) ) # Sigma    
    return   ( s.*((1.0.-t).^reverse(nn)) .* (t.^nn) )*p         
end

# function BezierFactorial2(p,n)
#     t   = collect(range(0,1.0,n))
#     N   = size(p,1)
#     nn  = collect(0:1:N-1)'
#     f   = factorial.(nn)
#     s   = factorial(N-1)./( f.*reverse(f) ) # Sigma
#     rnn = reverse(nn)
#     P   = zeros(n,3) # Curve points points
#     for q ∈ 1:1:n # Loop over all points
#         P[q,:] =  ( s.*((1.0-t[q]).^rnn) .* (t[q].^nn) )*p 
#     end
#     return P
# end

function BezierForm(p,n)
    t = collect(range(0,1.0,n))    
    
    c = 3.0*(p[2,:]-p[1,:])'    
    b = 3.0*(p[3,:]-p[2,:])'-c
    a = (p[4,:]-p[1,:])'.-c.-b

    P = a.*t.^3 .+ b.*t.^2 .+ c.*t .+ p[1,:]'.*ones(size(t))
    return P
end

# function BezierSet(p,n)

#     N = size(p,1) # Number of control points
#     i = 1    
#     while j<
#         j =


# end

# function BezierForm(p,n)
#     t = collect(range(0,1.0,n))
#     t²=t.^2
#     t³=t.^3
#     T1 = 1.0 .-3.0*t  .+3.0*t²      .-t³
#     T2 =       3.0*t  .-6.0*t²  .+3.0*t³
#     T3 =                3.0*t²  .-3.0*t³
#     T4 =                              t³
#     P = [T1 T2 T3 T4]*p
#     return P
# end

# function BezierFormN(p,n)
#     N = 
#     return P
# end


function BezierSpline(p,n)
    
    t   = collect(range(0,1.0,n)) # Parameter space vector
    
    # Characteristic matrix
    C = [ 1.0  0.0  0.0  0.0; 
         -3.0  3.0  0.0  0.0;
          3.0 -6.0  3.0  0.0;
         -1.0  3.0 -3.0  1.0 ] 

    return [ones(size(t)) t t.^2 t.^3]*C*p # Curve points
end

function CardinalSpline(p,n,s)
    
    t   = collect(range(0,1.0,n)) # Parameter space vector
    
    # Characteristic matrix
    C = [   0.0   1.0    0.0        0.0; 
             -s   0.0    s          0.0;
          2.0*s   s-3.0  3.0-2.0*s   -s;
             -s   2.0-s  s-2.0        s ] 

    return [ones(size(t)) t t.^2 t.^3]*C*p # Curve points
end

function CatmullRomSpline(p,n)
    
    t   = collect(range(0,1.0,n)) # Parameter space vector
    
    # Characteristic matrix
    C = [   0.0   1.0   0.0  0.0; 
           -0.5   0.0   0.5  0.0;
            1.0  -2.5   2.0 -0.5;
           -0.5   1.5  -1.5  0.5] 

    return [ones(size(t)) t t.^2 t.^3]*C*p # Curve points
end

function BSpline(p,n)
    
    t   = collect(range(0,1.0,n)) # Parameter space vector
    
    # Characteristic matrix
    C = [   1.0   4.0   1.0  0.0; 
           -3.0   0.0   3.0  0.0;
            3.0  -6.0   3.0  0.0;
           -1.0   3.0  -3.0  1.0]./6.0; 

    return [ones(size(t)) t t.^2 t.^3]*C*p # Curve points
end

P1 = BezierFactorial(p,n)
P2 = BezierForm(p,n)

P3 = CardinalSpline(p,n,0.5)
P4 = CatmullRomSpline(p,n)
P5 = BSpline(p,n)

function plotFunc(P,c)
    scatter!(ax, P[:,1], P[:,2], P[:,3],markersize=15,color=c)
    lines!(ax, P[:,1], P[:,2], P[:,3],linewidth=3,color=c)
end

# Visualize 
fig1 = Figure()

ax = Axis3(fig1[1, 1], aspect = :data)

plotFunc(p,:black)
plotFunc(P1,:blue)
plotFunc(P2,:red)
plotFunc(P3,:yellow)
plotFunc(P4,:green)
plotFunc(P5,:cyan)

fig1

