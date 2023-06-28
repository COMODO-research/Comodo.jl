# Call required packages

# using GLMakie # For visualisation
using GeometryBasics # For point and mesh format
using LinearAlgebra # For things like dot and cross products

# Function to create the surface geometry data for an icosahedron
function icosahedron(r=1.0)

    ϕ=(1.0+sqrt(5.0))/2.0 # Golden ratio

    # Create vertices
    s = r/sqrt(ϕ + 2.0)
    t = ϕ*s

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

    ϕ = (1.0+sqrt(5.0))/2.0 # Golden ratio
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
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,8)
    s = r/sqrt(3.0)
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

    # 4/3*b/a
    # V[1 ] = GeometryBasics.Point{3, Float64}(  -0.5, -b/6.0, -0.25*a/b)
    # V[2 ] = GeometryBasics.Point{3, Float64}(   0.5, -b/6.0, -0.25*a/b)    
    # V[3 ] = GeometryBasics.Point{3, Float64}(   0.0,    0.0,  0.75*a/b)
    # V[4 ] = GeometryBasics.Point{3, Float64}(   0.0,  b/3.0, -0.25*a/b)  

    # Create faces
    F = Vector{TriangleFace{Int64}}(undef,4)
    F[1 ] = TriangleFace{Int64}(1,2,3)
    F[2 ] = TriangleFace{Int64}(4,2,1)
    F[3 ] = TriangleFace{Int64}(4,3,2)
    F[4 ] = TriangleFace{Int64}(4,1,3)
    
    return GeometryBasics.Mesh(V,F)
end

function toGeometryBasicsMesh(VM,FM)

    # Loop over vertex matrix and convert to GeometryBasics vector of Points
    V=Vector{GeometryBasics.Point{3, Float64}}(undef,size(VM,1))
    for q ∈ 1:1:size(VM,1)
        V[q] = GeometryBasics.Point{3, Float64}(VM[q,:])
    end

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
    #else # Other mesh type

    end

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