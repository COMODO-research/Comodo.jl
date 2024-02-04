# Call required packages

using GeometryBasics # For point and mesh format
using LinearAlgebra # For things like dot and cross products
using DataStructures # For unique_dict
using Interpolations # E.g. for resampling curves
using Statistics # For: mean
using GLMakie # For slidercontrol
using SparseArrays # For meshconnectivity

function gibbondir()
    joinpath(@__DIR__, "..")
end

function slidercontrol(hSlider,ax)    
    on(events(ax).keyboardbutton) do event
        if event.action == Keyboard.press || event.action == Keyboard.repeat # Pressed or held for instance
            if event.key == Keyboard.up                                  
                sliderValue = hSlider.value.val              
                if sliderValue ==hSlider.range.val[end]                
                    set_close_to!(hSlider, hSlider.range.val[1])
                else                
                    set_close_to!(hSlider, hSlider.value.val+1)
                end            
            elseif event.key == Keyboard.down
                sliderValue = hSlider.value.val            
                if sliderValue ==hSlider.range.val[1]
                    set_close_to!(hSlider, hSlider.range.val[end])
                else                
                    set_close_to!(hSlider, hSlider.value.val-1)
                end                        
            end
            if event.key == Keyboard.right                                              
                set_close_to!(hSlider, hSlider.value.val+1)            
            elseif event.key == Keyboard.left            
                set_close_to!(hSlider, hSlider.value.val-1)
            end
        end
    end    
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

function gridpoints(X,Y=X,Z=X)    
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

function interp_biharmonic_spline(x,y,xi; extrapolate_method="linear",pad_data="linear")
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
            yi[.!L] = interp_biharmonic(xx,yy,xi[.!L]) # Use biharmonic interpolation for points within range
        else # Nothing to extrapolate
            yi = interp_biharmonic(xx,yy,xi)
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
            yi[L] = interp_biharmonic(xx,yy,xi[L]) # Use biharmonic interpolation for points within range
        else # Nothing to extrapolate
            yi = interp_biharmonic(xx,yy,xi)
        end
    elseif extrapolate_method=="biharmonic"
        # Allow extrapolation as per the biharmonic function
        yi = interp_biharmonic(xx,yy,xi) 
    end

    return yi
end

function interp_biharmonic(x,y,xi)
    # Distances from all points in X to all points in X
    Dxx = dist(x,x)

    # Determine weights for interpolation
    g = (Dxx.^2) .* (log.(Dxx).-1.0)   # Green's function.
    g[1:size(g,1)+1:length(g)] .= 0.0 # Fix values along diagonal
    g[isnan.(g)] .= 0.0 # Replace NaN entries by zeros 
    W = g \ y # Weights  

    D = dist(xi,x) # Distance between points in X and XI

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

function dist(V1,V2)
    D = Matrix{Float64}(undef,length(V1),length(V2))   
    for i ∈ eachindex(V1)
        for j ∈ eachindex(V2)          
            D[i,j] = norm((V1[i].-V2[j]))       
        end
    end
    return D
end

function mindist(V1,V2; getIndex=false, skipSelf = false )
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
function meshedges(M; unique_only=false)
    if isa(M,GeometryBasics.Mesh)
        S = faces(M) #Get simplices
    else
        S = M        
    end    
    n = length(S) #Number of simplices
    
    # Get number of nodes per simplex, use first for now (better to get this from some property)
    s1 = S[1] # First simplex
    m = length(s1) #Number of nodes per simplex
    
    E = GeometryBasics.LineFace{Int}[]    
    for j1 ∈ 1:1:m # Loop over each node/point for the current simplex           
        if j1<m
            j2 = j1+1
        else
            j2 = 1
        end            
        for s ∈ S # Loop over each simplex        
            push!(E,(s[j1],s[j2]))            
        end 
    end

    # Remove doubles e.g. 1-2 seen as same as 2-1
    if unique_only
        E = gunique(E; sort_entries=true);
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

function togeometrybasics_faces(FM)
    # Loop over face matrix and convert to GeometryBasics vector of Faces (e.g. QuadFace, or TriangleFace)
    n = length(FM)
    m = length(FM[1])
    if m == 3 # Triangles
        F = [TriangleFace{Int64}(FM[q]) for q ∈ eachindex(FM)]        
    elseif m ==4 # Quads
        F = [QuadFace{Int64}(FM[q]) for q ∈ eachindex(FM)]        
    else # Other mesh type        
        F = [NgonFace{Int64}(FM[q]) for q ∈ eachindex(FM)]        
    end
    return F
end

function togeometrybasics_faces(FM::Matrix{Int64})
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

function togeometrybasics_points(VM)
    # Loop over vertex matrix and convert to GeometryBasics vector of Points
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,length(VM))
    for q ∈ 1:1:size(VM,1)
        V[q] = GeometryBasics.Point{3, Float64}(VM[q])
    end
    return V
end

function togeometrybasics_points(VM::Matrix{Float64})
    # Loop over vertex matrix and convert to GeometryBasics vector of Points
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,size(VM,1))
    for q ∈ 1:1:size(VM,1)
        V[q] = GeometryBasics.Point{3, Float64}(VM[q,:])
    end
    return V
end

function togeometrybasics_points(VM::Vector{Vector{Int64}})
    # Loop over vertex matrix and convert to GeometryBasics vector of Points
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,length(VM))
    for q ∈ 1:1:size(VM,1)
        V[q] = GeometryBasics.Point{3, Float64}(VM[q])
    end
    return V
end

function togeometrybasics_mesh(VM,FM)

    V=togeometrybasics_points(VM)
    F=togeometrybasics_faces(FM)

    return GeometryBasics.Mesh(V,F)

end

function vertexnormal(M::GeometryBasics.Mesh) 
    F=faces(M) # Get faces
    V=coordinates(M) # Get vertices
    return vertexnormal(F,V)
end

function vertexnormal(F,V) 
    NF = facenormal(F,V)
    return vertexnormal(F,V)
end

function facenormal(M::GeometryBasics.Mesh) 
    F=faces(M) # Get faces
    V=coordinates(M) # Get vertices
    return facenormal(F,V)
end

# Computes mesh face normals
function facenormal(F,V)
    # Computes the normal vectors for the input surface geometry defined by the vertices V and the faces F
    N = Vector{GeometryBasics.Vec{3, Float64}}(undef,length(F)) # Allocate array for normal vectors
    n =  length(F[1]) # Number of nodes per face    
    for q ∈ eachindex(F) # Loop over all faces
        c  = cross(V[F[q][n]],V[F[q][1]]) # Initialise as cross product of last and first vertex position vector
        for qe ∈ 1:1:n-1 # Loop from first to end-1            
            c  += cross(V[F[q][qe]],V[F[q][qe+1]]) # Add next edge contribution          
        end
        N[q] = c./norm(c) # Normalise vector length and add
    end
    return N
end

function facenormal(M::GeometryBasics.Mesh) 
    F=faces(M) # Get faces
    V=coordinates(M) # Get vertices
    return facenormal(F,V)
end

# Computes mesh face normals
function facenormal(F,V)
    # Computes the normal vectors for the input surface geometry defined by the vertices V and the faces F
    N = Vector{GeometryBasics.Vec{3, Float64}}(undef,length(F)) # Allocate array for normal vectors
    n =  length(F[1]) # Number of nodes per face    
    for q ∈ eachindex(F) # Loop over all faces
        c  = cross(V[F[q][n]],V[F[q][1]]) # Initialise as cross product of last and first vertex position vector
        for qe ∈ 1:1:n-1 # Loop from first to end-1            
            c  += cross(V[F[q][qe]],V[F[q][qe+1]]) # Add next edge contribution          
        end
        N[q] = c./norm(c) # Normalise vector length and add
    end
    return N
end

function vertexnormal(M::GeometryBasics.Mesh) 
    NF = facenormal(M)
    V = coordinates(M)
    return simplex2vertexdata(faces(M),NF,V)
end

function vertexnormal(F,V) 
    NF = facenormal(F,V)
    return simplex2vertexdata(F,NF,V)
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
            indReverse[i] = j 
        else
            indReverse[i] = d[X[i]]            
        end        
    end
    return collect(keys(d)), indUnique, indReverse #, c
end

function unique_simplices(F,V=missing)
    if ismissing(V)
        n = maximum(reduce(vcat,F)) 
    else
        n = length(V)
    end

    virtualFaceIndices = sub2ind(n.*ones(Int64,length(F[1])),sort.(F))    
    ~, ind1, ind2 = unique_dict(virtualFaceIndices) 

    return F[ind1], ind1, ind2
end

function subtri(F,V,n; method = "linear")
    
    if n==0
        return F,V
    elseif n==1
        E = meshedges(F)
        Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        
        Fm1 = [TriangleFace{Int64}(a.+length(V)) for a ∈ eachrow(reshape(indReverse,length(F),length(F[1])))] 
        Fm2 = Vector{TriangleFace{Int64}}(undef,length(Fm1))
        Fm3 = Vector{TriangleFace{Int64}}(undef,length(Fm1))
        Fm4 = Vector{TriangleFace{Int64}}(undef,length(Fm1))        
        for i ∈ eachindex(F)                        
            Fm2[i] = TriangleFace{Int64}([Fm1[i][1], Fm1[i][3], F[i][1]])
            Fm3[i] = TriangleFace{Int64}([Fm1[i][2], Fm1[i][1], F[i][2]])
            Fm4[i] = TriangleFace{Int64}([Fm1[i][3], Fm1[i][2], F[i][3]])
        end

        # Create combined face set
        Fn = [Fm1; Fm2; Fm3; Fm4]        

        con_E2F = con_edge_face(F,Eu,indReverse)

        # Create new vertices depending on method
        if method == "linear" # Simple linear splitting
            # Create complete point set
            Vn = [V; simplexcenter(Eu,V)]  # Old and new mid-edge points          
        elseif method == "loop" #Loop subdivision 
    
            # New mid-edge like vertices
            Vm = Vector{GeometryBasics.Point{3, Float64}}(undef,length(Eu)) 
            for q ∈ eachindex(Eu) # For each edge index                        
                F_touch = F[con_E2F[q]] # Faces sharing current edge, mostly 2 but 1 for a boundary edge
                indVerticesTouch = Vector{Int64}() 
                for f ∈ F_touch        
                    b = f.!=Eu[q][1] .&& f.!=Eu[q][2]      
                    if any(b)  
                        append!(indVerticesTouch,f[b])           
                    end
                end        
                Vm[q]=3/8 .*(V[Eu[q][1]] .+ V[Eu[q][2]])  .+ 1/8 .* (V[indVerticesTouch[1]] .+ V[indVerticesTouch[2]])
            end
    
            # Modified vertices for original vertices
            Vv = Vector{GeometryBasics.Point{3, Float64}}(undef,length(V))
            for q ∈ eachindex(V)            
                B_vert_face = [any(f.==q) for f in F]
                F_touch = F[B_vert_face] # Faces mostly 2 but 1 for a boundary edge
                indVerticesTouch = Vector{Int64}()
                for f ∈ F_touch                
                    indTouch = f[f.!=q]        
                    for i ∈ indTouch 
                        if i ∉ indVerticesTouch 
                            push!(indVerticesTouch,i)
                        end
                    end
                end
                N = length(indVerticesTouch)                
                v_sum = sum(V[indVerticesTouch],dims=1)[1]                
                β = 1/N * (5/8-(3/8 +1/4*cos((2*π)/N))^2)        
                Vv[q] = (1-N*β) .* V[q] .+ β*v_sum                
            end    
            # Create complete point set
            Vn = [Vv;Vm] # Updated orignals and new "mid-edge-ish" points
        else
            error("""Incorrect metod. Use: "linear" or "loop" """)
        end

        return Fn,Vn    

    elseif n>1
        for _ =1:1:n
            F,V = subtri(F,V,1; method=method)
        end
        return F,V
    end
end

function subquad(F,V,n; method="linear")
    if n==0
        return F,V
    elseif n==1

        # Get edges
        E = meshedges(F) # Non-unique edges

        Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        
        con_F2E = con_face_edge(F,Eu,indReverse)
        con_E2F = con_edge_face(F,Eu,indReverse)
        con_V2E = con_vertex_edge(Eu,V)
        con_V2F = con_vertex_face(F,V)

        # Define vertices
        if method =="linear"
            Ve = simplexcenter(Eu,V) # Mid edge points
            Vf = simplexcenter(F,V)  # Mid face points
            Vn = [V;Ve;Vf] # Joined point set
        elseif method =="Catmull-Clark"
            # Mid face points
            Vf = simplexcenter(F,V)  
            Ve_mid = simplexcenter(Eu,V) # Mid edge points

            # Edge points 
            Ve = Vector{GeometryBasics.Point{3, Float64}}(undef,length(Eu)) # Initialize edge points
            for q ∈ eachindex(Eu)                         
                Ve[q] = (mean(Vf[con_E2F[q]],dims=1)[1] .+ Ve_mid[q])./2.0
            end

            # Vertex points 
            Vv = Vector{GeometryBasics.Point{3, Float64}}(undef,length(V)) # Initialize vertex points
            for q ∈ eachindex(V) # Loop over all vertices
                indF = con_V2F[q]
                indE = con_V2E[q]
                N = length(indF) # Number of faces (or edges) touching this vertex                    
                Vv[q] = (mean(Vf[indF],dims=1)[1] .+ 2.0.*mean(Ve_mid[indE],dims=1)[1] .+ (N-3.0).*V[q])./N
            end

            Vn = [Vv;Ve;Vf] # Joined point set
        end

        # Define faces
        Fn = Vector{QuadFace{Int64}}(undef,length(F)*4)        
        nv = length(V)
        ne = length(Eu)
        for q ∈ eachindex(F)
            i = 1 + (q-1)*4
            for ii = 0:1:3
                Fn[i+ii] = QuadFace{Int64}([F[q][ii+1], con_F2E[q][ii+1]+nv, q+nv+ne, con_F2E[q][1+mod(3+ii,4)]+nv])                
            end            
        end
        return Fn,Vn
    elseif n>1
        for _ =1:1:n
            F,V = subquad(F,V,1;method=method)
        end

        return F,V
    end
end

# Create geodesic dome
function geosphere(n,r)
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

function hexbox(boxDim,boxEl)    
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

mutable struct ConnectivitySet
    edge_vertex::Vector{LineFace{Int64}}
    edge_face::Vector{Vector{Int64}}
    edge_edge::Vector{Vector{Int64}}
    face_vertex#::Vector{Vector{Int64}}
    face_edge::Vector{Vector{Int64}}        
    face_face::Vector{Vector{Int64}}
    vertex_edge::Vector{Vector{Int64}}
    vertex_face::Vector{Vector{Int64}}        
    vertex_vertex::Vector{Vector{Int64}}
    ConnectivitySet(E_uni, con_E2F, con_E2E, F,  con_F2E, con_F2F, con_V2E, con_V2F, con_V2V) = new(E_uni, con_E2F, con_E2E, F,  con_F2E, con_F2F, con_V2E, con_V2F, con_V2V) 
end

function con_face_edge(F,E_uni=missing,indReverse=missing)
    if ismissing(E_uni) | ismissing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end
    return [Vector{Int64}(a) for a ∈ eachrow(reshape(indReverse,length(F),length(F[1])))] # [indReverse[[1,2,3].+ (i-1)*3] for i ∈ eachindex(F)]
end

function con_edge_face(F,E_uni=missing,indReverse=missing)
    if ismissing(E_uni)| ismissing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end
    indF = repeat(1:1:length(F); outer=length(F[1]))
    A = sparse(indReverse,indF,indF)
    con_E2F = Vector{Vector{Int64}}(undef,length(E_uni))
    for q ∈ eachindex(E_uni)
        _,v = findnz(A[q,:])
        con_E2F[q] = v
    end    
    return con_E2F
end

function con_face_face(F,E_uni=missing,indReverse=missing,con_E2F=missing,con_F2E=missing)
    if ismissing(E_uni)| ismissing(indReverse)
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
    end
    
    if ismissing(con_E2F) 
        con_E2F = con_edge_face(F,E_uni)
    end

    if ismissing(con_F2E)
        con_F2E = con_face_edge(F,E_uni,indReverse)                 
    end

    con_F2F = Vector{Vector{Int64}}(undef,length(F))
    for q ∈ eachindex(F)
        s = sparsevec(reduce(vcat,con_E2F[con_F2E[q]]),1)
        s[q] = 0 # Null and thus remove the "self" entry
        con_F2F[q] = findall(!iszero, s)
    end

    return con_F2F
end

function con_vertex_face(F,V=missing)
    if ismissing(V)
        n = maximum(reduce(vcat,F))
    else
        n = length(V)
    end
    indF = repeat(1:1:length(F); inner=length(F[1]))
    B = sparse(reduce(vcat,F),indF,indF)
    con_V2F = Vector{Vector{Int64}}(undef,n)
    for q ∈ 1:1:n
        _,v = findnz(B[q,:])
        con_V2F[q] = v
    end
    return con_V2F
end
function con_vertex_edge(E_uni,V=missing)
    if ismissing(V)
        n = maximum(reduce(vcat,E_uni))
    else 
        n = length(V)
    end
    indE = repeat(1:1:length(E_uni); inner=length(E_uni[1]))
    C = sparse(reduce(vcat,E_uni),indE,indE)
    con_V2E = Vector{Vector{Int64}}(undef,n)
    for q ∈ 1:1:n
        _,v = findnz(C[q,:])
        con_V2E[q] = v
    end
    return con_V2E
end

function con_edge_edge(E_uni,con_V2E=missing)
    if ismissing(con_V2E)
        con_V2E = con_vertex_edge(E_uni) 
    end    
    con_E2E = Vector{Vector{Int64}}(undef,length(E_uni))
    for q ∈ eachindex(E_uni)
        s = sparsevec(reduce(vcat,con_V2E[E_uni[q]]),1)
        s[q] = 0 # Null and thus remove the "self" entry
        con_E2E[q] = findall(!iszero, s)
    end
    return con_E2E
end

function con_vertex_vertex(E_uni,V=missing,con_V2E=missing)
    if ismissing(V)
        n = maximum(reduce(vcat,E_uni))
    else 
        n = length(V)
    end
    if ismissing(con_V2E)
        con_V2E = con_vertex_edge(E_uni,V)
    end

    con_V2V = Vector{Vector{Int64}}(undef,n)
    for q ∈ 1:1:n
        s = sparsevec(reduce(vcat,E_uni[con_V2E[q]]),1)
        s[q] = 0 # Null and thus remove the "self" entry
        con_V2V[q] = findall(!iszero, s)
    end
    return con_V2V
end

function meshconnectivity(F,V) #conType = ["ev","ef","ef","fv","fe","ff","ve","vf","vv"]

    # EDGE-VERTEX connectivity
    E = meshedges(F)
    E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)

    # FACE-EDGE connectivity
    con_F2E = con_face_edge(F,E_uni,indReverse)    

    # EDGE-FACE connectivity
    con_E2F = con_edge_face(F,E_uni)

    # FACE-FACE connectivity
    con_F2F = con_face_face(F,E_uni,indReverse,con_E2F,con_F2E)

    # VERTEX-FACE connectivity
    con_V2F = con_vertex_face(F,V)

    # VERTEX-EDGE connectivity
    con_V2E = con_vertex_edge(E_uni,V)

    # EDGE-EDGE connectivity
    con_E2E = con_edge_edge(E_uni,con_V2E)

    # VERTEX-VERTEX connectivity
    con_V2V = con_vertex_vertex(E_uni,V,con_V2E)

    return ConnectivitySet(E_uni, con_E2F, con_E2E, F,  con_F2E, con_F2F, con_V2E, con_V2F, con_V2V)
end 

function mergevertices(F,V; roundVertices = true, numDigitsMerge=missing)

    m = length(V)
    if roundVertices
        if ismissing(numDigitsMerge)
            E = meshedges(F)
            d = [sqrt( sum((V[e[1]] .- V[e[2]]).^2) ) for e ∈ E]
            pointSpacing = mean(d)
            m = round(Int64,log10(pointSpacing))
            numDigitsMerge = 6-m
        end

        # Create rounded coordinates to help obtain unique set
        VR = [round.(v,digits = numDigitsMerge) for v ∈ V]

        # Get unique indices and reverse for rounded vertices
        _,ind1,ind2 = gunique(VR; return_index=true, return_inverse=true,sort_entries=false)
        V = V[ind1] # The unique node set
    else
        V,ind1,ind2 = gunique(V; return_index=true, return_inverse=true,sort_entries=false)
    end

    if length(V) != m # If the length has changed
        # Correct indices for faces
        for q ∈ eachindex(F)
            F[q] = ind2[F[q]]
        end
    end

    return F,V,ind1,ind2
end

function smoothmesh_laplacian(F,V,con_V2V=missing; n=1, λ=0.5)
    """
    This function implements Weighted Laplacian mesh smoothing. At each 
    iteration, this method replaces each point by the weighted sum of the 
    Laplacian mean for the point and the point itself. The weighting is 
    controlled by the parameter λ which is in the range (0,1). If λ=0 no 
    smoothing occurs. If λ=0 pure Laplacian mean based smoothing occurs. For 
    intermediate values a linear blending between the two occurs.     
    """
    if λ!=1
        # Compute vertex-vertex connectivity i.e. "Laplacian umbrellas" if missing
        if ismissing(con_V2V)
            E = meshedges(F)
            E_uni,_,_ = unique_simplices(E)
            con_V2V = con_vertex_vertex(E_uni)
        end
        for _ = 1:1:n
            Vs = deepcopy(V)
            for q ∈ eachindex(V)
                # Linear blend between original and pure Laplacian
                Vs[q] = (1.0-λ).*Vs[q] .+ λ*mean(V[con_V2V[q]])
            end
            V = Vs
        end
    end
    return V
end

function smoothmesh_hc(F,V, con_V2V=missing; n=1, α=0.1, β=0.5, tolDist=missing)
    """
    This function implements HC (Humphrey's Classes) smoothing. This method uses
    Laplacian like smoothing but aims to compensate for the shrinkage/swelling 
    seen with pure Laplacian smoothing. 
    
    REF: 
    Vollmer et al. Improved Laplacian Smoothing of Noisy Surface Meshes, 1999
    https://doi.org/10.1111/1467-8659.00334
    """
    # Compute vertex-vertex connectivity i.e. "Laplacian umbrellas" if missing
    if ismissing(con_V2V)
        E = meshedges(F)
        E_uni,_,_ = unique_simplices(E)
        con_V2V = con_vertex_vertex(E_uni)
    end
    P = deepcopy(V) # Copy original input points
    B = deepcopy(V) # Initialise B
    c = 0
    while c<n       
        Q = deepcopy(P) # Reset Q as P for this iteration
        for i ∈ eachindex(V)
            P[i] = mean(Q[con_V2V[i]]) # Laplacian 
            # Compute different vector between P and a point between original 
            # point and Q (which is P before laplacian)
            B[i] = P[i] .- (α.*V[i] .+ (1-α).*Q[i])
        end
        d = 0.0        
        for i ∈ eachindex(V)      
            # Push points back based on blending between pure difference vector
            # B and the Laplacian mean of these      
            P[i] = P[i] .- (β.*B[i] .+ (1-β).* mean(B[con_V2V[i]]))
        end   
        c+=1 
        if !ismissing(tolDist) # Include tolerance based termination
            d = 0.0
            for i ∈ eachindex(V)
                d+=sqrt(sum((P[i].-Q[i]).^2)) # Sum of distances
            end
            if d<tolDist # Sum of distance smaller than tolerance?
                break
            end            
        end
    end
    return P
end

function quadplate(plateDim,plateElem)
    num_x = plateElem[1]+1
    num_y = plateElem[2]+1
    V = Vector{GeometryBasics.Point{3, Float64}}()
    for y = range(-plateDim[2]/2,plateDim[2]/2,num_y)
        for x = range(-plateDim[1]/2,plateDim[1]/2,num_x)
            push!(V,GeometryBasics.Point{3, Float64}(x,y,0.0))
        end
    end

    F = Vector{QuadFace{Int64}}()
    ij2ind(i,j) = i + ((j-1)*num_x) # function to convert subscript to linear indices

    for i = 1:1:plateElem[1]
        for j = 1:1:plateElem[2]        
            push!(F,QuadFace{Int64}([ij2ind(i,j),ij2ind(i+1,j),ij2ind(i+1,j+1),ij2ind(i,j+1)]))
        end
    end
    return F, V
end

function quadsphere(r,n)
    M = platonicsolid(2,r)
    F = faces(M)
    V = coordinates(M)
    if n > 0
        for q ∈ 1:1:n
            F,V = subquad(F,V,1)
            V = r .* (V ./ norm.(V))
        end
    end
    return F,V
end

function simplex2vertexdata(F,DF,V=missing)
    con_V2F = con_vertex_face(F,V)
    DV = (typeof(DF))(undef,length(con_V2F))
    T = eltype(DV)
    for q ∈ eachindex(DV)
        DV[q] = mean(T,DF[con_V2F[q]])
    end
    return DV
end

function vertex2simplexdata(F,DV)
    T = eltype(DV) # Element type of data in DV
    DF =  (typeof(DV))(undef,length(F)) # Allocate data for F
    for q ∈ eachindex(F)
        DF[q] = mean(T,DV[F[q]]) # The mean of the vertex data for each entry in F
    end
    return DF
end

function simplexcenter(F,V)
    return vertex2simplexdata(F,V)
end

function normalizevector(A)
    return A./norm.(A)
end

function circlepoints(r,n)
    return [GeometryBasics.Point{3, Float64}(r*cos(t),r*sin(t),0) for t ∈ range(0,2*π-(2*π)/n,n)]
end


function loftlinear(V1,V2;num_loft=2,close_loop=true)

    num_loop = length(V1)
    T = eltype(V1)
    # Linearly blending points from first to last
    V = Vector{T}()
    for q ∈ range(0,num_loft,num_loft)
        λ = q/num_loft
        Vn = (1.0-λ).*V1 .+ λ.* V2 
        append!(V,Vn)
    end

    # Build faces
    F = Vector{QuadFace{Int64}}()
    ij2ind(i,j) = i + ((j-1)*num_loop) # function to convert subscript to linear indices

    for i = 1:1:num_loop-1
        for j = 1:1:num_loft-1        
            push!(F,QuadFace{Int64}([ij2ind(i,j),ij2ind(i+1,j),ij2ind(i+1,j+1),ij2ind(i,j+1)]))
        end
    end

    # Add faces to close over shape if requested
    if close_loop
        for q ∈ 1:1:num_loft-1
            push!(F,QuadFace{Int64}([ij2ind(1,q),ij2ind(1,q+1),ij2ind(num_loop,q+1),ij2ind(num_loop,q)]))
        end
    end

    return F, V
end 

function normalplot(ax,M; type_flag="face", color=:black,linewidth=2)
    F = faces(M)
    V = coordinates(M)
    E = meshedges(F)
    d = mean([norm(V[e[1]]-V[e[2]]) for e ∈ E])
    if type_flag == "face"        
        N = facenormal(F,V)
        V = simplexcenter(F,V)        
    elseif type_flag == "vertex"
        N = vertexnormal(F,V)          
    else
        error(""" Incorrect type_flag, use "face" or "vertex" """)
    end 
    # hp = arrows!(ax,V,N.*d,color=color,quality=6)    
    E = [GeometryBasics.LineFace{Int}(i,i+length(V)) for i ∈ 1:1:length(V)]
    V = append!(V,V.+N.*d)
    hp = wireframe!(ax,GeometryBasics.Mesh(V,E),linewidth=linewidth, transparency=false, color=color)
    return hp 
end



