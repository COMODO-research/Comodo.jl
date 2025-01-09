using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Random

Random.seed!(1)

function randangle(siz)
    A = Matrix{Float64}(undef,siz)    
    for i in eachindex(A)
        A[i] = rand()*pi*rand((-1,1))
    end
    return A
end


# Define grid
size_grid = (25,25) # grid size
sampleFactor = 35 # Pixel sample factor wrt grid
pixelSize = 1/sampleFactor # Pixel size assuming grid has unit steps

# Create grid vectors 
A = randangle(size_grid) # Random angles

# for w = 0:1:19
# A .+= (2*pi/20)

# A = pi/2 .* ones(size_grid)
# for i in eachindex(A)
#     A[i] *= rand((-1,1))
#     # println(i)
# end
# A[7]*=-1

Ux = cos.(A) # Unit vector x-component
Uy = sin.(A) # Unit vector y-component

# Initialise image
size_image = (size_grid .- 1) .* sampleFactor # image size
M = Matrix{Float64}(undef,size_image) # Start as undef Float64

# Pre-compute grid cell quantities
xy = range(0+pixelSize/2,1-pixelSize/2,sampleFactor) # x or y coordinates within a grid cell

X = [x for x in 0:sampleFactor:(size_grid[2]-1)*sampleFactor, j in 0:sampleFactor:(size_grid[1]-1)*sampleFactor]'
Y = [y for i in 0:sampleFactor:(size_grid[2]-1)*sampleFactor, y in 0:sampleFactor:(size_grid[1]-1)*sampleFactor]'

# ig = ceil(Int,i./sampleFactor) # Grid row index
                # ip = mod1(i,sampleFactor) # Cell row index

                # jg = ceil(Int,j./sampleFactor) # Grid column index
                # jp = mod1(j,sampleFactor) # Cell column index

# ceil.(collect(1:1:25)./sampleFactor)

function fade_perlin(a,b,t)
    # q = 3.0*t^3 - 2.0*t^3 # Smoothstep
    q = 0.5 .- 0.5.*cos.(t*pi) # Cosine
    # q = t # Linear 

    # Fade function q = 6t⁵-15t⁴+10t³    
    # q = t^3 * (t * (6.0 * t - 15.0) + 10.0) # Perlin fade function
    
    return (1.0-q)*a + q*b
end

xc = [0,1,1,0]
yc = [0,0,1,1]
for ip in 1:sampleFactor # For each cell row
    for jp in 1:sampleFactor # For each cell column
        for ig in 1:size_grid[1]-1 # For each grid row    
            for jg in 1:size_grid[2]-1 # For each grid column
                i = (ig-1)*sampleFactor + ip # Pixel row index
                j = (jg-1)*sampleFactor + jp # Pixel column index
                
                # Current pixel cell coordinates
                px = xy[jp]
                py = xy[ip]
                        
                # Offset vector components
                xc1 = px # -xc[1] Offset vector 1 x
                xc2 = px-xc[2] # Offset vector 2 x
                xc3 = px-xc[3] # Offset vector 3 x
                xc4 = px-xc[4] # Offset vector 4 x

                yc1 = py # -yc[2] Offset vector 1 y                
                yc2 = py-yc[2] # Offset vector 2 y               
                yc3 = py-yc[3] # Offset vector 3 y                
                yc4 = py-yc[4] # Offset vector 4 y

                u1x = Ux[ig  ,jg]
                u2x = Ux[ig  ,jg+1]
                u3x = Ux[ig+1,jg+1]
                u4x = Ux[ig+1,jg]

                u1y = Uy[ig  ,jg]
                u2y = Uy[ig  ,jg+1]
                u3y = Uy[ig+1,jg+1]
                u4y = Uy[ig+1,jg]

                d1 = xc1.*u1x + yc1.*u1y
                d2 = xc2.*u2x + yc2.*u2y
                d3 = xc3.*u3x + yc3.*u3y
                d4 = xc4.*u4x + yc4.*u4y
        
                # Interpolation 
                # Fade function 6t⁵-15t⁴+10t³
                d12 = fade_perlin(d1,d2,px)
                d34 = fade_perlin(d4,d3,px)
                d = fade_perlin(d12,d34,py)

                M[i,j] = d
                
            end
        end
    end
end

fig = Figure(size=(1500,1500))
# ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Edge angles")
ax1 = Axis(fig[1, 1], aspect = DataAspect(), title = "Perlin noise",limits=(-sampleFactor,size_image[2]+sampleFactor,-sampleFactor,size_image[1]+sampleFactor) )
hm = image!(ax1, M',interpolate=false,colormap = Makie.Reverse(:Spectral),colorrange=(-0.5,0.5)) #

arrows!(ax1,reduce(vcat,X), reduce(vcat,Y), reduce(vcat,Ux), reduce(vcat,Uy), arrowsize = 10, lengthscale = sampleFactor/3, color = :black, linewidth=1)

Colorbar(fig[1, 2], hm)

fig

# save("/home/kevin/DATA/Julia/Comodo.jl/assets/img/perlin_noise_c"*string(w)*".jpg",fig,px_per_unit = 1)
# save("/home/kevin/DATA/Julia/Comodo.jl/assets/img/perlin_noise_large.jpg",fig,px_per_unit = 2)

# end