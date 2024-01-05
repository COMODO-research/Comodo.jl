# Call required packages

using GeometryBasics # For point and mesh format
using LinearAlgebra # For things like dot and cross products
using DataStructures # For unique_dict
using Interpolations # E.g. for resampling curves

function gibbonDir()
    joinpath(@__DIR__, "..")
end

function elements2indices(F)
    n = length(F) # Number of elements
    m = length(F[1]) # Number of nodes per element
    ind = Vector{Int64}(undef,n*m)
    for i ∈ eachindex(F)
        ii = 1 + (i-1)*m
        ind[ii:ii+m-1] = F[i] #convert(Vector{Int64},F[i])
    end
    return unique(ind)
end

function rot3(α=0.0,β=0.0,γ=0.5*pi)
    # Creates a rotation tensor (matrix) based on the input Euler angles α, β, and γ

    # Rotation around x-axis
    Qx=[ 1.0     0.0     0.0
         0.0     cos(α) -sin(α)
         0.0     sin(α)  cos(α)]

    # Rotation around y-axis
    Qy=[ cos(β)  0.0     sin(β)
         0.0     1.0     0.0
        -sin(β)  0.0     cos(β)]

    # Rotation around z-axis
    Qz=[ cos(γ) -sin(γ)  0.0
         sin(γ)  cos(γ)  0.0
         0.0     0.0     1.0]
        
    return Qx*Qy*Qz # Return combined rotation tensor
end

function gridPoints(X,Y=X,Z=X)    
    # Creates a 3D grid of GeometryBasics type 3D points based on the iterable input defining the ranges in the x, y, and z direction
    Vg = Vector{GeometryBasics.Point{3, Float64}}(undef,length(X)*length(Y)*length(Z)) # Allocate grid point set
    q=0 # Initialize index 
    for x ∈ X # For all x entries
        for y ∈ Y # For all y entries
            for z ∈ Z # For all z entries
                q+=1 # Increment index
                Vg[q]=GeometryBasics.Point{3, Float64}([x,y,z]) # Add the current point
            end
        end
    end
    return Vg # Return the grid point set
end

function interp_biharmonicSpline(x,y,xi; extrapolate_method="linear",pad_data="linear")
    # This function uses biharmonic spline interpolation. The input is assumed to represent ordered data representing a curve
    # 
    # Reference:  David T. Sandwell, Biharmonic spline interpolation of 
    # GEOS-3 and SEASAT altimeter data, Geophysical Research Letters, 2,
    # 139-142, 1987.  

    # Pad data if needed
    xx = collect(x)
    yy = collect(y)
    if pad_data=="linear"
        # Linearly extended ends are added 
        dx1 = x[1]-x[2]
        dy1 = y[1]-y[2]    
        dx2 = x[end]-x[end-1]
        dy2 = y[end]-y[end-1]
           
        pushfirst!(xx,x[1]+dx1)
        pushfirst!(xx,x[1]+2.0*dx1) 
        push!(xx,x[end]+dx2)
        push!(xx,x[end]+2.0*dx2)
    
        pushfirst!(yy,y[1]+dy1)
        pushfirst!(yy,y[1]+2.0*dy1) 
        push!(yy,y[end]+dy2)
        push!(yy,y[end]+2.0*dy2)
    elseif pad_data=="constant"
        # The start and end are copied twice 
        dx1 = x[1]-x[2]       
        dx2 = x[end]-x[end-1]
            
        pushfirst!(xx,x[1]+dx1)
        pushfirst!(xx,x[1]+2.0*dx1) 
        push!(xx,x[end]+dx2)
        push!(xx,x[end]+2.0*dx2)
    
        pushfirst!(yy,y[1])
        pushfirst!(yy,y[1]) 
        push!(yy,y[end])
        push!(yy,y[end])
    elseif  pad_data=="none"
        # No padding
    else
        error(""" Invalid pad_data method provided, valued options are: "linear", "constant", and "none" """)
    end

    # Change behaviour depending on extrapolation method
    if extrapolate_method=="linear"
        # Simple data based linear extrapolation 
        L = [xii<x[1] || xii>x[end] for xii in xi] # Boolean for points to extrapolate for
        if any(L) # If any points outside of the range were encountered
            yi = Vector{Float64}(undef,length(xi)) # Initialise yi
            yi[L] = lerp(xx,yy,xi[L]) # Linearly extrapolate outside of the input range
            yi[.!L] = interp_biharmonic_ND(xx,yy,xi[.!L]) # Use biharmonic interpolation for points within range
        else # Nothing to extrapolate
            yi = interp_biharmonic_ND(xx,yy,xi)
        end
    elseif extrapolate_method=="constant"
        # Simple constant extrapolation (last/first value is repeated indefinately)
        Ls = xi.<x[1] # Boolean for points preceeding the start
        Ll = xi.>x[end] # Boolean for points after the end
        if any(Ls .|| Ll) # If any points outside of the range were encountered
            yi = Vector{Float64}(undef,length(xi)) # Initialise yi
            if any(Ls) # If any before start were found 
                yi[Ls] .= yy[1] # Just copy the start for these
            end
            if any(Ll) # If any after the end were found
                yi[Ll] .= yy[end] # Just copy the end for these
            end
            L = .!Ls .&& .!Ll # Boolean for points to interpolate using biharmonic interpolation
            yi[L] = interp_biharmonic_ND(xx,yy,xi[L]) # Use biharmonic interpolation for points within range
        else # Nothing to extrapolate
            yi = interp_biharmonic_ND(xx,yy,xi)
        end
    elseif extrapolate_method=="biharmonic"
        # Allow extrapolation as per the biharmonic function
        yi = interp_biharmonic_ND(xx,yy,xi) 
    end

    return yi
end

function interp_biharmonic_ND(x,y,xi)
    # Distances from all points in X to all points in X
    Dxx = distND(x,x)

    # Determine weights for interpolation
    g = (Dxx.^2) .* (log.(Dxx).-1.0)   # Green's function.
    g[1:size(g,1)+1:length(g)] .= 0.0 # Fix values along diagonal
    g[isnan.(g)] .= 0.0 # Replace NaN entries by zeros 
    W = g \ y # Weights  

    D = distND(xi,x) # Distance between points in X and XI

    G = (D.^2).*(log.(D).-1.0) # Green's function.
    G[isnan.(G)] .= 0.0 # Replace NaN entries by zeros
    return G * W
end

function lerp(x,y,xi)    
    if length(xi)>1
        yi=zeros(eltype(xi),length(xi))
        for q ∈ eachindex(xi)
            yi[q] = lerp_(x,y,xi[q])
        end
    else
        yi = lerp_(x,y,xi)
    end 
    return yi
end

function lerp_(x,y,xi)
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
    
    w = abs(x[j]-x[i])
    t = (xi-x[i])/w
    yi = (1.0-t)*y[i] + t*y[j]
    return yi
end

function distND(V1,V2)
    D = Matrix{Float64}(undef,length(V1),length(V2))   
    for i ∈ eachindex(V1)
        for j ∈ eachindex(V2)          
            D[i,j] = sqrt(sum((V1[i].-V2[j]).^2))       
        end
    end
    return D
end

function minDist(V1,V2; getIndex=false, skipSelf = false )
    D = Vector{Float64}(undef,length(V1))
    d = Vector{Float64}(undef,length(V2))
    if getIndex
        I = Vector{Int64}(undef,length(V1))
    end
    for i ∈ eachindex(V1)
        for j ∈ eachindex(V2)
            if skipSelf && i==j
                d[j] = Inf
            else
                d[j] = sqrt(sum((V1[i].-V2[j]).^2)) 
            end       
        end
        if getIndex
            D[i], I[i] = findmin(d)
        else
            D[i] = minimum(d)
        end
    end
    if getIndex
        return D, I
    else
        return D
    end
end

function unique_dict_index(X; sort_entries=false)
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    T = eltype(X)
    d = Dict{T,Nothing}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indUnique = Vector{Int64}() 
    for i ∈ eachindex(X)        
        if sort_entries && length(X[1])>1
            x = sort(X[i])
        else
            x = X[i]
        end
        if !haskey(d, x)
            d[x] = nothing
            push!(xUni, X[i]) #Grow unique set
            push!(indUnique, i) #Grow unique indices
        end
    end
    return xUni, indUnique
end

function unique_dict_index_inverse(X; sort_entries=false)
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indUnique = Vector{Int64}() 
    indInverse = Vector{Int64}(undef,length(X)) 
    j=0
    for i ∈ eachindex(X)         
        if sort_entries && length(X[1])>1
            x = sort(X[i])        
        else
            x = X[i]
        end
        if !haskey(d, x)
            j+=1 # Increment counter
            d[x] = j # inverse index in dict
            push!(xUni, X[i]) #Grow unique set
            push!(indUnique, i) #Grow unique indices
            indInverse[i] = j  # Store inverse index          
        else
            indInverse[i]=d[x]
        end
    end
    return xUni, indUnique, indInverse
end

function unique_dict_index_count(X; sort_entries=false)
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indUnique = Vector{Int64}() 
    c =  Vector{Int64}() 

    j=0
    for i ∈ eachindex(X)      
        if sort_entries && length(X[1])>1
            x = sort(X[i])        
        else
            x = X[i]
        end  
        if !haskey(d, x)
            j+=1 # Increment counter
            d[x] = j # inverse index in dict
            push!(xUni, X[i]) #Grow unique set
            push!(indUnique, i) #Grow unique indices
            push!(c,1)
        else
            c[d[x]] += 1
        end
    end
    return xUni, indUnique, c
end


function unique_dict_index_inverse_count(X; sort_entries=false)
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indUnique = Vector{Int64}() 
    indInverse = Vector{Int64}(undef,length(X)) 
    c =  Vector{Int64}() 

    j=0
    for i ∈ eachindex(X)      
        if sort_entries && length(X[1])>1
            x = sort(X[i])        
        else
            x = X[i]
        end  
        if !haskey(d, x)
            j+=1 # Increment counter
            d[x] = j # inverse index in dict
            push!(xUni, X[i]) #Grow unique set
            push!(indUnique, i) #Grow unique indices
            indInverse[i] = j  # Store inverse index
            push!(c,1)
        else
            indInverse[i]=d[x]
            c[d[x]] += 1
        end
    end
    return xUni, indUnique, indInverse,c
end

function unique_dict_count(X; sort_entries=false)
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    c = Vector{Int64}()
    j = 0
    for i ∈ eachindex(X)      
        if sort_entries && length(X[1])>1
            x = sort(X[i])        
        else
            x = X[i]
        end  
        if !haskey(d, x)
            j += 1 # Index in unique array
            d[x] = j # Store count in dict
            push!(xUni, X[i]) # Grow unique set
            push!(c,1) # Grow counting array
        else
            c[d[x]] += 1
        end
    end
    return xUni, c
end

function unique_dict_inverse(X; sort_entries=false)
    # Here a normal Dict is used to keep track of unique elements. Normal dicts do not maintain element insertion order. 
    # Hence the unique indices need to seperately be tracked. 
    T = eltype(X)
    d = Dict{T,Int64}() # Use dict to keep track of used values
    xUni = Vector{T}()
    indInverse = Vector{Int64}(undef,length(X)) 

    j=0
    for i ∈ eachindex(X)      
        if sort_entries && length(X[1])>1
            x = sort(X[i])        
        else
            x = X[i]
        end  
        if !haskey(d, x)
            j+=1 # Increment counter
            d[x] = j # inverse index in dict
            push!(xUni, X[i]) #Grow unique set
            indInverse[i] = j  # Store inverse index
        else
            indInverse[i]=d[x]
        end
    end
    return xUni, indInverse
end 


function gunique(X; return_unique=true, return_index=false, return_inverse=false, return_counts=false, sort_entries=false)
    # Use required unique function 
    if return_unique==true && return_index==false && return_inverse==false && return_counts==false
        # UNIQUE
        if sort_entries && length(X[1])>1
            return unique(sort.(X))
        else
            return unique(X)
        end
    elseif return_unique==true && return_index==true && return_inverse==false && return_counts==false
        # UNIQUE, INDICES
        return unique_dict_index(X; sort_entries=sort_entries)
    elseif return_unique==true && return_index==false && return_inverse==false && return_counts==true
        # UNIQUE, COUNTS
        return unique_dict_count(X; sort_entries=sort_entries)
    elseif return_unique==true && return_index==false && return_inverse==true && return_counts==false
        # UNIQUE, INVERSE
        return unique_dict_inverse(X; sort_entries=sort_entries)
    elseif return_unique==true && return_index==true && return_inverse==true && return_counts==false
        # UNIQUE, INDICES, INVERSE
        return unique_dict_index_inverse(X; sort_entries=sort_entries)
    elseif return_unique==true && return_index==true && return_inverse==false && return_counts==true
        # UNIQUE, INDICES, COUNTS
        return unique_dict_index_count(X; sort_entries=sort_entries)
    elseif return_unique==true && return_index==true && return_inverse==true && return_counts==true
        # UNIQUE, INDICES, INVERSE, COUNTS
        return unique_dict_index_inverse_count(X; sort_entries=sort_entries)
    end
end

"""
    ind2sub(siz,ind)
Converts the linear indices in `ind`, for a matrix/array with size `siz`, to the equivalent subscript indices.  
"""

function ind2sub(siz,ind)
    
    numDim = length(siz)
    k = cumprod([siz[i] for i ∈ eachindex(siz)],dims=1)
    m = prod(siz)
    
    if any(ind.>m) || any(ind.<1)
        error("Encountered index value of out of valid range 1:"*string(m))
    end

    if length(ind)>1
        A = [ind2sub_(ind_i,numDim,k) for ind_i ∈ ind]
    else
        A = ind2sub_(ind,numDim,k)      
    end

    return A
end

# ind2sub helper function to parse just a single linear index and produce a single subscript index set 
function ind2sub_(ind,numDim,k)

    a = Vector{Int64}(undef,numDim) # Initialise a

    for q ∈ numDim:-1:1   # For all dimensions     
        if q==1 # First 1st dimension
            a[1] = rem(ind-1,k[1]) + 1        
        else       
            p = rem(ind-1,k[q-1]) + 1 # "previous"
            a[q] = (ind - p)./k[q-1] + 1 # Current        
            ind = p # Set indices as "previous"          
        end            
    end
    
    return a
end

"""
    sub2ind(siz,A)
Converts the subscript indices in `A`, for a matrix/array with size `siz`, to the equivalent linear indices.  
"""

function sub2ind(siz,A)
    
    numDim = length(siz)
    k = cumprod([siz[i] for i ∈ eachindex(siz)],dims=1)

    ind = Vector{Int64}(undef,length(A))
    for i ∈ eachindex(A)        
        ind[i] = round(Int64,A[i][1]);
        for q=2:1:numDim
            ind[i] += (A[i][q].-1).*k[q-1]
        end        
    end

    return ind
end

function sub2ind(siz,A::Vector{Int64})
    
    numDim = length(siz)
    k = cumprod([siz[i] for i ∈ eachindex(siz)],dims=1)
    
    ind = A[1];
    for q=2:1:numDim
        ind += (A[q].-1).*k[q-1]
    end        

    return ind
end

# Function to obtain the edges of a face based mesh
function meshEdges(M)
    if isa(M,GeometryBasics.Mesh)
        S = faces(M) #Get simplices
    else
        S = M        
    end    
    n=length(S) #Number of simplices
    
    # Get number of nodes per simplex, use first for now (better to get this from some property)
    s1 = S[1]
    m=length(s1) #Number of nodes per simplex
    
    # E = Vector{Vector{Int64}}(undef,n*m)  
    E = [Vector{Int64}(undef,2) for _ in 1:n*m]
    i=1
    for q ∈ eachindex(S)
        s = S[q]        
        for j1 ∈ 1:1:m            
            if j1<m
                j2=j1+1
            else
                j2=1
            end            
            E[i]=[s[j1],s[j2]]
            i+=1
        end 
    end
    return E
end

# Function to create the surface geometry data for an icosahedron
function icosahedron(r=1.0)

    ϕ = Base.MathConstants.golden # (1.0+sqrt(5.0))/2.0, Golden ratio
    s = r/sqrt(ϕ + 2.0)
    t = ϕ*s

    # Create vertices
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,12)
    V[1 ] = GeometryBasics.Point{3, Float64}( 0.0,   -s,  -t)
    V[2 ] = GeometryBasics.Point{3, Float64}( 0.0,   -s,   t)
    V[3 ] = GeometryBasics.Point{3, Float64}( 0.0,    s,   t)
    V[4 ] = GeometryBasics.Point{3, Float64}( 0.0,    s,  -t)
    V[5 ] = GeometryBasics.Point{3, Float64}(  -s,   -t, 0.0)
    V[6 ] = GeometryBasics.Point{3, Float64}(  -s,    t, 0.0)
    V[7 ] = GeometryBasics.Point{3, Float64}(   s,    t, 0.0)
    V[8 ] = GeometryBasics.Point{3, Float64}(   s,   -t, 0.0)
    V[9 ] = GeometryBasics.Point{3, Float64}(  -t,  0.0,  -s)
    V[10] = GeometryBasics.Point{3, Float64}(   t,  0.0,  -s)
    V[11] = GeometryBasics.Point{3, Float64}(   t,  0.0,   s)
    V[12] = GeometryBasics.Point{3, Float64}(  -t,  0.0,   s)

    # Create faces
    F = Vector{TriangleFace{Int64}}(undef,20)
    F[1 ] = TriangleFace{Int64}(9,4,1)
    F[2 ] = TriangleFace{Int64}(1,5,9)
    F[3 ] = TriangleFace{Int64}(1,8,5)
    F[4 ] = TriangleFace{Int64}(10,8,1)
    F[5 ] = TriangleFace{Int64}(4,10,1)
    F[6 ] = TriangleFace{Int64}(5,2,12)
    F[7 ] = TriangleFace{Int64}(12,2,3)
    F[8 ] = TriangleFace{Int64}(12,3,6)
    F[9 ] = TriangleFace{Int64}(12,6,9)
    F[10] = TriangleFace{Int64}(12,9,5)
    F[11] = TriangleFace{Int64}(10,7,11)
    F[12] = TriangleFace{Int64}(8,10,11)
    F[13] = TriangleFace{Int64}(2,8,11)
    F[14] = TriangleFace{Int64}(3,2,11)
    F[15] = TriangleFace{Int64}(7,3,11)
    F[16] = TriangleFace{Int64}(2,5,8)
    F[17] = TriangleFace{Int64}(10,4,7)
    F[18] = TriangleFace{Int64}(7,6,3)
    F[19] = TriangleFace{Int64}(6,7,4)
    F[20] = TriangleFace{Int64}(6,4,9)
    
    return GeometryBasics.Mesh(V,F)
end

# Function to create the surface geometry data for a octahedron
function octahedron(r=1.0)

    s = r/sqrt(2.0)

    # Create vertices    
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,6)
    V[1 ] = GeometryBasics.Point{3, Float64}(   -s,    -s,  0.0)
    V[2 ] = GeometryBasics.Point{3, Float64}(    s,    -s,  0.0)
    V[3 ] = GeometryBasics.Point{3, Float64}(    s,     s,  0.0)
    V[4 ] = GeometryBasics.Point{3, Float64}(   -s,     s,  0.0)
    V[5 ] = GeometryBasics.Point{3, Float64}(  0.0,   0.0,   -r)
    V[6 ] = GeometryBasics.Point{3, Float64}(  0.0,   0.0,    r)
    
    # Create faces
    F = Vector{TriangleFace{Int64}}(undef,8)
    F[1 ] = TriangleFace{Int64}(5,2,1)
    F[2 ] = TriangleFace{Int64}(5,3,2)
    F[3 ] = TriangleFace{Int64}(5,4,3)
    F[4 ] = TriangleFace{Int64}(5,1,4)
    F[5 ] = TriangleFace{Int64}(6,1,2)
    F[6 ] = TriangleFace{Int64}(6,2,3)
    F[7 ] = TriangleFace{Int64}(6,3,4)
    F[8 ] = TriangleFace{Int64}(6,4,1)
    
    return GeometryBasics.Mesh(V,F)
end



# Function to create the surface geometry data for a dodecahedron
function dodecahedron(r=1.0)

    ϕ = Base.MathConstants.golden # (1.0+sqrt(5.0))/2.0, Golden ratio
    s = r/sqrt(3.0)
    t = ϕ*s    
    w = (ϕ-1.0)*s

    # Create vertices    
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,20)
    V[1 ] = GeometryBasics.Point{3, Float64}(   s,   s,   s)
    V[2 ] = GeometryBasics.Point{3, Float64}(   w, 0.0,   t)
    V[3 ] = GeometryBasics.Point{3, Float64}(  -t,  -w, 0.0)
    V[4 ] = GeometryBasics.Point{3, Float64}(   t,   w, 0.0)
    V[5 ] = GeometryBasics.Point{3, Float64}(  -s,   s,  -s)
    V[6 ] = GeometryBasics.Point{3, Float64}( 0.0,  -t,  -w)
    V[7 ] = GeometryBasics.Point{3, Float64}(  -t,   w, 0.0)
    V[8 ] = GeometryBasics.Point{3, Float64}(   s,  -s,   s)
    V[9 ] = GeometryBasics.Point{3, Float64}(  -s,   s,   s)
    V[10] = GeometryBasics.Point{3, Float64}(  -s,  -s,   s)
    V[11] = GeometryBasics.Point{3, Float64}(   s,  -s,  -s)
    V[12] = GeometryBasics.Point{3, Float64}(   w, 0.0,  -t)
    V[13] = GeometryBasics.Point{3, Float64}(  -s,  -s,  -s)
    V[14] = GeometryBasics.Point{3, Float64}( 0.0,  -t,   w)
    V[15] = GeometryBasics.Point{3, Float64}( 0.0,   t,  -w)
    V[16] = GeometryBasics.Point{3, Float64}(  -w, 0.0,   t)
    V[17] = GeometryBasics.Point{3, Float64}(   t,  -w, 0.0)
    V[18] = GeometryBasics.Point{3, Float64}(  -w, 0.0,  -t)
    V[19] = GeometryBasics.Point{3, Float64}(   s,   s,  -s)
    V[20] = GeometryBasics.Point{3, Float64}( 0.0,   t,   w)

    # Create faces
    F = Vector{NgonFace{5,Int64}}(undef,12)
    F[1 ] = NgonFace{5,Int64}(20, 9,16, 2, 1)
    F[2 ] = NgonFace{5,Int64}( 2,16,10,14, 8)
    F[3 ] = NgonFace{5,Int64}(16, 9, 7, 3,10)
    F[4 ] = NgonFace{5,Int64}( 7, 9,20,15, 5)
    F[5 ] = NgonFace{5,Int64}(18,13, 3, 7, 5)
    F[6 ] = NgonFace{5,Int64}( 3,13, 6,14,10)
    F[7 ] = NgonFace{5,Int64}( 6,13,18,12,11)
    F[8 ] = NgonFace{5,Int64}( 6,11,17, 8,14)
    F[9 ] = NgonFace{5,Int64}(11,12,19, 4,17)
    F[10] = NgonFace{5,Int64}( 1, 2, 8,17, 4)
    F[11] = NgonFace{5,Int64}( 1, 4,19,15,20)
    F[12] = NgonFace{5,Int64}(12,18, 5,15,19)
    
    return GeometryBasics.Mesh(V,F)
end

# Function to create the surface geometry data for a cube
function cube(r=1.0)

    # Create vertices       
    s = r/sqrt(3.0)

    V=Vector{GeometryBasics.Point{3, Float64}}(undef,8)
    V[1 ] = GeometryBasics.Point{3, Float64}( -s, -s, -s)
    V[2 ] = GeometryBasics.Point{3, Float64}( -s,  s, -s)
    V[3 ] = GeometryBasics.Point{3, Float64}(  s,  s, -s)
    V[4 ] = GeometryBasics.Point{3, Float64}(  s, -s, -s)
    V[5 ] = GeometryBasics.Point{3, Float64}( -s, -s,  s)
    V[6 ] = GeometryBasics.Point{3, Float64}( -s,  s,  s)
    V[7 ] = GeometryBasics.Point{3, Float64}(  s,  s,  s)
    V[8 ] = GeometryBasics.Point{3, Float64}(  s, -s,  s)
        
    # Create faces
    F = Vector{QuadFace{Int64}}(undef,6)
    F[1 ] = QuadFace{Int64}(1,2,3,4)
    F[2 ] = QuadFace{Int64}(8,7,6,5)
    F[3 ] = QuadFace{Int64}(5,6,2,1)
    F[4 ] = QuadFace{Int64}(6,7,3,2)    
    F[5 ] = QuadFace{Int64}(7,8,4,3)    
    F[6 ] = QuadFace{Int64}(8,5,1,4)    

    return GeometryBasics.Mesh(V,F)
end

# Function to create the surface geometry data for a tetrahedron
function tetrahedron(r=1.0)

    # Create vertices    
    a = r*sqrt(2.0)/sqrt(3.0)
    b = -r*sqrt(2.0)/3.0
    c = -r/3.0       

    V=Vector{GeometryBasics.Point{3, Float64}}(undef,4)
    V[1 ] = GeometryBasics.Point{3, Float64}(   -a,      b,   c)
    V[2 ] = GeometryBasics.Point{3, Float64}(    a,      b,   c)    
    V[3 ] = GeometryBasics.Point{3, Float64}(  0.0,    0.0,   r)
    V[4 ] = GeometryBasics.Point{3, Float64}(  0.0, -2.0*b,   c)  

    # Create faces
    F = Vector{TriangleFace{Int64}}(undef,4)
    F[1 ] = TriangleFace{Int64}(1,2,3)
    F[2 ] = TriangleFace{Int64}(4,2,1)
    F[3 ] = TriangleFace{Int64}(4,3,2)
    F[4 ] = TriangleFace{Int64}(4,1,3)
    
    return GeometryBasics.Mesh(V,F)
end

function toGeometryBasicsSimplices(FM)
    # Loop over face matrix and convert to GeometryBasics vector of Faces (e.g. QuadFace, or TriangleFace)
    if size(FM,2)==3 # Triangles
        F = Vector{TriangleFace{Int64}}(undef,size(FM,1))
        for q ∈ 1:1:size(FM,1)            
            F[q] = TriangleFace{Int64}(FM[q,:])
        end
    elseif size(FM,2)==4 # Quads
        F = Vector{QuadFace{Int64}}(undef,size(FM,1))
        for q ∈ 1:1:size(FM,1)            
            F[q] = QuadFace{Int64}(FM[q,:])
        end
    else # Other mesh type
        m=size(FM,2)
        F = Vector{m,NgonFace{Int64}}(undef,size(FM,1))        
        for q ∈ 1:1:size(FM,1)            
            F[q] = NgonFace{m,Int64}(FM[q,:])
        end
    end
    return F
end

function toGeometryBasicsPoints(VM)
    # Loop over vertex matrix and convert to GeometryBasics vector of Points
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,size(VM,1))
    for q ∈ 1:1:size(VM,1)
        V[q] = GeometryBasics.Point{3, Float64}(VM[q,:])
    end
    return V
end

function toGeometryBasicsPoints(VM::Vector{Vector{Float64}})
    # Loop over vertex matrix and convert to GeometryBasics vector of Points
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,size(VM,1))
    for q ∈ 1:1:size(VM,1)
        V[q] = GeometryBasics.Point{3, Float64}(VM[q])
    end
    return V
end

function toGeometryBasicsPoints(VM::Vector{Vector{Int64}})
    # Loop over vertex matrix and convert to GeometryBasics vector of Points
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,size(VM,1))
    for q ∈ 1:1:size(VM,1)
        V[q] = GeometryBasics.Point{3, Float64}(VM[q])
    end
    return V
end

function toGeometryBasicsMesh(VM,FM)

    V=toGeometryBasicsPoints(VM)
    F=toGeometryBasicsSimplices(FM)

    return GeometryBasics.Mesh(V,F)

end


function meshnormal(M::GeometryBasics.Mesh) 
    F=faces(M) # Get faces
    V=coordinates(M) # Get vertices

    # Computes the normal vectors for the input surface geometry defined by the vertices V and the faces F
    N = Vector{GeometryBasics.Point{3, Float64}}(undef,length(F)) # Allocate array for normal vectors
    VN = Vector{GeometryBasics.Point{3, Float64}}(undef,length(F))  # Allocate mid-face coordinates
    n =  length(F[1]) # Number of nodes per face    
    for q ∈ eachindex(F) # Loop over all faces
        c  = cross(V[F[q][n]],V[F[q][1]]) # Initialise as cross product of last and first vertex position vector
        vn = V[F[q][n]]./n # Initialise mean face coordinate computation (hence division by n)
        for qe ∈ 1:1:n-1 # Loop from first to end-1            
            c  += cross(V[F[q][qe]],V[F[q][qe+1]]) # Add next edge contribution          
            vn += V[F[q][qe]]./n # Add next vertex contribution
        end
        a = sqrt(sum(c.^2)) # Compute vector magnitude
        N[q] = c./a # Normalise vector length and add
        VN[q] = vn # Add mean face coordinate        
    end
    return N, VN
end

# Computes mesh face normals
function meshnormal(V::Matrix{Float64},F::Matrix{Int64})
    # Computes the normal vectors for the input surface geometry defined by the vertices V and the faces F
    N = zeros(Float64,size(F,1),3) # Allocate array for normal vectors
    VN = zeros(Float64,size(F,1),3) # Allocate mid-face coordinates
    n = size(F,2) # Number of nodes per face
    for q ∈ 1:1:size(F,1) # Loop over all faces
        c  = cross(V[F[q,n],:],V[F[q,1],:]) # Initialise as cross product of last and first vertex position vector
        vn = V[F[q,n],:]./n # Initialise mean face coordinate computation (hence division by n)
        for qe ∈ 1:1:n-1 # Loop from first to end-1            
            c  += cross(V[F[q,qe],:],V[F[q,qe+1],:]) # Add next edge contribution          
            vn += V[F[q,qe],:]./n # Add next vertex contribution
        end
        a = sqrt(sum(c.^2)) # Compute vector magnitude
        N[q,:] = c./a # Normalise vector length and add
        VN[q,:] = vn # Add mean face coordinate
    end
    return N, VN
end

# Creates mesh data for a platonic solid of choice
function platonicsolid(n,r=1.0)

    if n==1
        M = tetrahedron(r)
    elseif n==2
        M = cube(r)
    elseif n==3
        M = octahedron(r)
    elseif n==4 
        M = icosahedron(r)
    elseif n==5
        M = dodecahedron(r)
    end

    return M
end

function unique_dict(X)    
    d = OrderedDict{eltype(X),Int64}() # Use dict to keep track of used values    
    indUnique = Vector{Int64}()
    indReverse = Vector{Int64}(undef,length(X))
    # c = Vector{Int64}()
    j=0
    for i ∈ eachindex(X)        
        if !haskey(d, X[i])                                          
            j+=1
            d[X[i]] = j # reverse index in dict            
            push!(indUnique, i)                     
            # push!(c, 1) 
            indReverse[i] = j 
        else
            indReverse[i] = d[X[i]]            
            # c[d[X[i]]] += 1
        end        
    end
    return collect(keys(d)), indUnique, indReverse #, c
end

function unique_simplices(F,V)        
    virtualFaceIndices = sub2ind(length(V).*ones(Int64,length(F[1])),sort.(F))    
    ~, ind1, ind2 = unique_dict(virtualFaceIndices) 

    return F[ind1], ind1, ind2
end

function midEdgePoints(E,V)
    Vm = Vector{GeometryBasics.Point{3, Float64}}(undef,length(E)) 
    for q ∈ eachindex(E)
        Vm[q]=0.5.*(V[E[q][1]].+V[E[q][2]])
    end
    return Vm
end

function subtri(F,V,n)
    
    if n==0
        return F,V
    elseif n==1
        E = meshEdges(F)
        numEdges = size(E,1)
        Eu, ~, ind2 = unique_simplices(E,V)  
        Vm = midEdgePoints(Eu,V)
        Fmm = reshape(ind2,3,Int64(numEdges/3))'.+length(V)
        Fm1 = toGeometryBasicsSimplices(Fmm)

        Fm2 = Vector{TriangleFace{Int64}}(undef,length(Fm1))
        Fm3 = Vector{TriangleFace{Int64}}(undef,length(Fm1))
        Fm4 = Vector{TriangleFace{Int64}}(undef,length(Fm1))        
       for i ∈ eachindex(F)                        
            Fm2[i] = TriangleFace{Int64}([Fm1[i][1], Fm1[i][3], F[i][1]])
            Fm3[i] = TriangleFace{Int64}([Fm1[i][2], Fm1[i][1], F[i][2]])
            Fm4[i] = TriangleFace{Int64}([Fm1[i][3], Fm1[i][2], F[i][3]])
        end
        # Join new point and face sets
        Vn = [V;Vm]
        Fn = [Fm1; Fm2; Fm3; Fm4]        
        return Fn,Vn

    elseif n>1
        for _ =1:1:n
            F,V = subtri(F,V,1)
        end
        return F,V
    end
end

# Create geodesic dome
function geoSphere(n,r)
    M = platonicsolid(4,r)
    V = coordinates(M)
    F = faces(M)
    for _ = 1:n
        F,V = subtri(F,V,1)
        for q ∈ eachindex(V)
            v = V[q]
            rn = sqrt(sum(v.^2))
            V[q] = v .* (r/rn)
        end
    end
    return F,V
end

function hexMeshBox(boxDim,boxEl)    
    boxNod = boxEl.+1 # Number of nodes in each direction
    numElements = prod(boxEl) # Total number of elements
    numNodes = prod(boxNod) # Total number of nodes

    # Create hexahedral element description 
    ind1 = 1:numElements
    ijk_shift = [[0,0,0], 
                [1,0,0],
                [1,1,0],
                [0,1,0],
                [0,0,1],
                [1,0,1],
                [1,1,1],
                [0,1,1]]
                
    E = [Vector{Int64}(undef,8) for _ in 1:numElements] # Allocate elements

    for q ∈ 1:numElements
        ijk_1 = ind2sub(boxEl,ind1[q])    
        ijk_2 = ijk_1 .+ ijk_shift[2]
        ijk_3 = ijk_1 .+ ijk_shift[3]
        ijk_4 = ijk_1 .+ ijk_shift[4]
        ijk_5 = ijk_1 .+ ijk_shift[5]
        ijk_6 = ijk_1 .+ ijk_shift[6]
        ijk_7 = ijk_1 .+ ijk_shift[7]
        ijk_8 = ijk_1 .+ ijk_shift[8]
        
        E[q] = [sub2ind(boxNod,ijk_1),sub2ind(boxNod,ijk_2),sub2ind(boxNod,ijk_3),sub2ind(boxNod,ijk_4),sub2ind(boxNod,ijk_5),sub2ind(boxNod,ijk_6),sub2ind(boxNod,ijk_7),sub2ind(boxNod,ijk_8)]
    end

    # Create vertices aka nodal coordinates
    indNodes = collect(Int64,1:numNodes)
    IJK_nodes = ind2sub(boxNod,indNodes)

    V = convert(Vector{GeometryBasics.Point{3, Float64}},IJK_nodes)
    for q ∈ eachindex(V)
        V[q]=(V[q].-[1,1,1]).*(boxDim./boxEl)
    end

    # Create face sets from elements
    F = Vector{QuadFace{Int64}}(undef,numElements*6) # Allocate faces
    CF_type = Vector{Int64}(undef,numElements*6) # Allocate face color/label data
    for q = eachindex(E) # Loop over all elements
        e = E[q] # The current element 
        qf = 1 + (q-1)*6 # Index mapping into face array 

        # Add the current element's faces
        F[qf  ] = QuadFace{Int64}(e[4],e[3],e[2],e[1]) # Top
        F[qf+1] = QuadFace{Int64}(e[5],e[6],e[7],e[8]) # Bottom
        F[qf+2] = QuadFace{Int64}(e[1],e[2],e[6],e[5]) # Side 1
        F[qf+3] = QuadFace{Int64}(e[3],e[4],e[8],e[7]) # Side 2
        F[qf+4] = QuadFace{Int64}(e[2],e[3],e[7],e[6]) # Front
        F[qf+5] = QuadFace{Int64}(e[4],e[1],e[5],e[8]) # Back

        # Add the current element's face color/labels
        CF_type[qf  ] = 1 # Top
        CF_type[qf+1] = 2 # Bottom
        CF_type[qf+2] = 3 # Side 1
        CF_type[qf+3] = 4 # Side 2
        CF_type[qf+4] = 5 # Front
        CF_type[qf+5] = 6 # Back    
    end

    F_uni,indUni,c_uni = gunique(F,return_index=true, return_inverse=false,return_counts=true,sort_entries=true)
    Lb = c_uni.==1
    Fb = F_uni[Lb]
    CF_type_uni = CF_type[indUni]
    CFb_type = CF_type_uni[Lb]

    return E,V,F,Fb,CFb_type
end