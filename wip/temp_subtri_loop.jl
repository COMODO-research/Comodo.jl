using Gibbon, GLMakie, GeometryBasics

function subtri_loop(F,V,n)
    
    if n==0
        return F,V
    elseif n==1
        
        E = meshEdges(F)
        numEdges = size(E,1)
        Eu, ~, ind2 = unique_simplices(E,V)  
        Fmm = reshape(ind2,3,Int64(numEdges/3))'
    
        Fm1 = toGeometryBasicsSimplices(Fmm.+length(V))
        Fm2 = Vector{TriangleFace{Int64}}(undef,length(Fm1))
        Fm3 = Vector{TriangleFace{Int64}}(undef,length(Fm1))
        Fm4 = Vector{TriangleFace{Int64}}(undef,length(Fm1))        
        for i ∈ eachindex(F)                        
            Fm2[i] = TriangleFace{Int64}([Fm1[i][1], Fm1[i][3], F[i][1]])
            Fm3[i] = TriangleFace{Int64}([Fm1[i][2], Fm1[i][1], F[i][2]])
            Fm4[i] = TriangleFace{Int64}([Fm1[i][3], Fm1[i][2], F[i][3]])
        end
    
        Vm = Vector{GeometryBasics.Point{3, Float64}}(undef,length(Eu)) 
        for q ∈ eachindex(Eu) # For each edge index                        
            F_touch = F[indIn(Eu[q],F)] # Faces sharing current edge, mostly 2 but 1 for a boundary edge
            indVerticesTouch = Vector{Int64}() 
            for f ∈ F_touch        
                b = f.!=Eu[q][1] .&& f.!=Eu[q][2]      
                if any(b)  
                    append!(indVerticesTouch,f[b])           
                end
            end
    
            Vm[q]=3/8 .*(V[Eu[q][1]] .+ V[Eu[q][2]])  .+ 1/8 .* (V[indVerticesTouch[1]] .+ V[indVerticesTouch[2]])
        end
    
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
    
        Vn = [Vv;Vm]
        Fn = [Fm1; Fm2; Fm3; Fm4]        
        return Fn,Vn    
    elseif n>1
        for _ =1:1:n
            F,V = subtri_loop(F,V,1)
        end
        return F,V
    end
end


r = 1.0 #radius
M = platonicsolid(4,r)
V = coordinates(M)
F = faces(M)

################# 

fig = Figure(size=(800,800))

stepRange = 0:1:4
hSlider = Slider(fig[2, 1], range = stepRange, startvalue = 0,linewidth=30)

titleString = lift(hSlider.value) do stepIndex
    "n = " * string(stepIndex) * " refinement steps"
end

Mn = lift(hSlider.value) do stepIndex
    Fn,Vn = subtri_loop(F,V,stepIndex)
    return GeometryBasics.Mesh(Vn,Fn)
end

ax = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString)
sliderControl(hSlider,ax)

hp1=wireframe!(ax,M,linewidth=3,color=:red, overdraw=false)
hp2=poly!(ax,Mn,strokewidth=2,color=:white, shading = FastShading)

# hp2=poly!(scene,M,strokewidth=1,color=:white, shading = FastShading)

fig