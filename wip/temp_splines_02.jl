using GLMakie
# using BenchmarkTools
# using Makie.GeometryBasics
# using ColorSchemes

n = 25

# Define raw data
x = -2:1:8 # Interval definition
y = 5.0*cos.(x.^2 ./ 9.0) 

xi = -2:0.5:8 #range(-2,8,n) # Interval definition




function central_diff(x,y)
    dx = diff(x)
    dy = diff(y)
    u = zeros(Float64,length(x))
    v = zeros(Float64,length(x))
    for q ∈ eachindex(x)
        if q==1
            u[1] = dx[1]
            v[1] = dy[1]
        elseif q==length(x)
            u[end] = dx[end]
            v[end] = dy[end]
        else
            u[q] = 0.5*(dx[q] + dx[q-1])
            v[q] = 0.5*(dy[q] + dy[q-1])
        end
    end
    return u,v
end

function v_spline(x,y)
    u = zeros(Float64,length(x))
    v = zeros(Float64,length(x))
    for q ∈ eachindex(x)
        if q==1
            u[1] = x[2]-x[1]
            v[1] = y[2]-y[1]
        elseif q==length(x)
            u[end] = x[end]-x[end-1]
            v[end] = y[end]-y[end-1]
        else
            u[q] = x[q+1] - x[q-1]
            v[q] = y[q+1] - y[q-1]
        end
    end
    return u,v
end

u,v = v_spline(x,y)

p = [0 0 0; 1 1 1; 2 -1 1; 3 0 0;] # Control points



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

function hermiteSpline(p,t)
        
    # Characteristic matrix
    C = [   1.0   0.0   0.0  0.0; 
            0.0   1.0   0.0  0.0;
           -3.0  -2.0   3.0 -1.0;
            2.0   1.0  -2.0  1.0]; 

    return [ones(size(t)) t t.^2 t.^3]*C*p # Curve points
end

function interpHermite(x,y,xi)
    vx,vy=v_spline(x,y) # Velocity vectors

    if length(xi)>1
        yi=zeros(eltype(xi),length(xi))
        xii=zeros(eltype(xi),length(xi))
        for q ∈ eachindex(xi)
            a,b = interpHermite_(x,y,vx,vy,xi[q])       
            yi[q] =a[1]
            xii[q]=b[1]        
            println("xi="*string(xi[q]))
            println("xt="*string(xii[q]))
        end
    else
        yi,xii = interpHermite_(x,y,vx,vy,xi)
    end 
    return yi,xii
end

function interpHermite_(x,y,vx,vy,xi)
    j = findfirst(x.>xi)
    if j==1 # Deal with extrapolation at the start
        i=1
        j=2
    elseif isnothing(j) # Deal with extrapolation at the end
        j = length(x)
        i = j-1 
    else
        i = j-1
    end
    
    w = x[j]-x[i]
    t = (xi-x[i])/w
    p = [x[i] y[i]; vx[i]/w vy[i]/w; x[j] y[j]; vx[j]/w vy[j]/w]
    p_i = hermiteSpline(p,t)

    println("i="*string(i))
    println("j="*string(j))
    println("t="*string(t))
    return p_i[:,2], p_i[:,1]
end

y_i,x_i = interpHermite(x,y,xi)
# P1 = BezierFactorial([x y],n)
# P2 = BezierForm(p,n)

# P3 = CardinalSpline(p,n,0.5)
# P4 = CatmullRomSpline(p,n)
# P5 = BSpline(p,n)

function plotFunc(x,y,c,s)
    scatter!(ax, x,y,markersize=s,color=c)
    lines!(ax, x,y,linewidth=3,color=c)
end

# Visualize 
fig1 = Figure()

ax = Axis(fig1[1, 1], aspect = DataAspect())

plotFunc(x,y,:black,25)
arrows!(ax,x, y, u, v,color=:red)
plotFunc(xi,y_i,:green,15)
plotFunc(x_i,y_i,:blue,10)

fig1

