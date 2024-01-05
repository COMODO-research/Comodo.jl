using GLMakie
using Interpolations

# Define raw data
x = range(-2,8,6) # Interval definition
y = 5.0*cos.(x.^2 ./ 9.0) 
n = 25
xi = range(-2,8,n) # Interval definition

function to_t(x,xi)
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

    return i,j,t,w 
end

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

function hermiteSpline(p,t)
        
    # Characteristic matrix
    C = [   1.0   0.0   0.0  0.0; 
            0.0   1.0   0.0  0.0;
           -3.0  -2.0   3.0 -1.0;
            2.0   1.0  -2.0  1.0]; 

    return [ones(size(t)) t t.^2 t.^3]*C*p # Curve points
end

function interp_cubic_hermite(x,y,xi)
    u,v = central_diff(x,y)
    yi = zeros(eltype(xi),length(xi))
    xii = zeros(eltype(xi),length(xi))

    for q ∈ eachindex(xi)
        i,j,t,w = to_t(x,xi[q])
        t3=t.^3
        t2=t.^2
        
        h00 = ( 2.0.*t3 .- 3.0.*t2           .+ 1.0)
        h10 = ( 1.0.*t3 .- 2.0.*t2 .+ 1.0.*t       )
        h01 = (-2.0.*t3 .+ 3.0.*t2                 )
        h11 = ( 1.0.*t3 .- 1.0.*t2                 )

        p1 = [x[i] y[i]]
        p2 = [x[j] y[j]]
        m1 = [u[i] v[i]]./w
        m2 = [u[j] v[j]]./w
        p_i = h00.*p1 .+ h10.*m1 .+ h01.*p2 + h11.*m2 
        yi[q]=p_i[2]
        xii[q]=p_i[1]
    end

    return yi,xii
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
    i,j,t,w = to_t(x,xi)
    p = [x[i] y[i]; vx[i]/w vy[i]/w; x[j] y[j]; vx[j]/w vy[j]/w]
    p_i = hermiteSpline(p,t)

    println("i="*string(i))
    println("j="*string(j))
    println("t="*string(t))
    return p_i[:,2], p_i[:,1]
end


u,v = v_spline(x,y)
uc,vc = central_diff(x,y)



y_i,x_i = interpHermite(x,y,xi)
y_i2,x_i2 = interp_cubic_hermite(x,y,xi)
itp_cubic = cubic_spline_interpolation(x, y, extrapolation_bc = Interpolations.Line())
y_i3 = itp_cubic(xi)

function plotFunc(x,y,c,s)
    scatter!(ax, x,y,markersize=s,color=c)
    lines!(ax, x,y,linewidth=3,color=c)
end

# Visualize 
fig1 = Figure()

ax = Axis(fig1[1, 1], aspect = DataAspect())

plotFunc(x,y,:black,25)
arrows!(ax,x, y, u, v,color=:red)
arrows!(ax,x, y, uc, vc,color=:blue)

plotFunc(xi,y_i,:green,15)
plotFunc(xi,y_i2,:blue,10)
plotFunc(xi,y_i3,:yellow,10)

fig1

