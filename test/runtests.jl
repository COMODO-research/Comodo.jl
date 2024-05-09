using Test, FileIO, Comodo, Comodo.GeometryBasics, Statistics, LinearAlgebra, GLMakie, Rotations

# ConnectivitySet

@testset "comododir" begin
    f = comododir()
    @test any(contains.(readdir(f),"src"))
    @test any(contains.(readdir(f),"assets"))
    @test any(contains.(readdir(f),"test"))
end

@testset "slidercontrol" verbose = true begin
    r = range(-2,2,10)
    startvalue = r[1]
    fig = Figure(size=(800,800))
    ax1 = Axis3(fig[1, 1])
    hSlider = Slider(fig[2, 1], range = r,startvalue = startvalue,linewidth=30)
    fig

    slidercontrol(hSlider,ax1)

    @test hSlider.selected_index[] == 1 # Test that slidercontrol did not alter index
end


@testset "slider2anim" verbose=true begin
    r = range(-2,2,10)
    startvalue = r[1]
    fig = Figure(size=(800,800))
    ax1 = Axis3(fig[1, 1])
    hSlider = Slider(fig[2, 1], range = r,startvalue = startvalue,linewidth=30)
    fig
    
    fileName = comododir()*"/assets/temp_anim1.mp4"
    slider2anim(fig,hSlider,fileName; backforth=true, duration=2)
    @test isfile(fileName) # File exists 
    rm(fileName) # Clean up

    fileName = comododir()*"/assets/temp_anim2.mp4"
    slider2anim(fig,hSlider,fileName; backforth=false, duration=pi)
    @test isfile(fileName) # File exists 
    rm(fileName) # Clean up

    fileName = comododir()*"/assets/temp_anim3.gif"
    slider2anim(fig,hSlider,fileName; backforth=false, duration=1.0)
    @test isfile(fileName) # File exists 
    rm(fileName) # Clean up
end


@testset "elements2indices" verbose = true begin
    @testset "Tri. faces" begin
        F = Vector{TriangleFace{Int64}}(undef, 3)
        F[1] = TriangleFace{Int64}(9, 4, 1)
        F[2] = TriangleFace{Int64}(1, 5, 9)
        F[3] = TriangleFace{Int64}(1, 8, 5)
        result = elements2indices(F)
        @test sort(result) == [1, 4, 5, 8, 9]
    end

    @testset "Quad. faces" begin
        F = Vector{QuadFace{Int64}}(undef, 6)
        F[1] = QuadFace{Int64}(1, 2, 3, 4)
        F[2] = QuadFace{Int64}(8, 7, 6, 5)
        F[3] = QuadFace{Int64}(5, 6, 2, 1)
        F[4] = QuadFace{Int64}(6, 7, 3, 2)
        F[5] = QuadFace{Int64}(7, 8, 4, 3)
        F[6] = QuadFace{Int64}(8, 5, 1, 4)
        result = elements2indices(F)
        @test sort(result) == [1, 2, 3, 4, 5, 6, 7, 8]
    end
end


@testset "gridpoints" verbose = true begin

    @testset "using ranges" begin
        # Define ranges of different types
        a = 1:2
        b = 1:1:2
        c = range(0,0.75,3)
        V = gridpoints(a,b,c)
        
        @test V == Point3{Float64}[[1.0, 1.0, 0.0], [2.0, 1.0, 0.0], 
        [1.0, 2.0, 0.0], [2.0, 2.0, 0.0], [1.0, 1.0, 0.375], 
        [2.0, 1.0, 0.375], [1.0, 2.0, 0.375], [2.0, 2.0, 0.375],
         [1.0, 1.0, 0.75], [2.0, 1.0, 0.75], [1.0, 2.0, 0.75], [2.0, 2.0, 0.75]]      
    end

    @testset "with 1 vector" begin
        a = Float64[1, 2, 3]

        expected = Point3{Float64}[[1.0, 1.0, 1.0], [2.0, 1.0, 1.0], 
        [3.0, 1.0, 1.0], [1.0, 2.0, 1.0], [2.0, 2.0, 1.0], [3.0, 2.0, 1.0], 
        [1.0, 3.0, 1.0], [2.0, 3.0, 1.0], [3.0, 3.0, 1.0], [1.0, 1.0, 2.0], 
        [2.0, 1.0, 2.0], [3.0, 1.0, 2.0], [1.0, 2.0, 2.0], [2.0, 2.0, 2.0], 
        [3.0, 2.0, 2.0], [1.0, 3.0, 2.0], [2.0, 3.0, 2.0], [3.0, 3.0, 2.0], 
        [1.0, 1.0, 3.0], [2.0, 1.0, 3.0], [3.0, 1.0, 3.0], [1.0, 2.0, 3.0], 
        [2.0, 2.0, 3.0], [3.0, 2.0, 3.0], [1.0, 3.0, 3.0], [2.0, 3.0, 3.0], 
        [3.0, 3.0, 3.0]]

        result = gridpoints(a)

        @test result == expected
    end

    @testset "with 2 vectors" begin
        a = Float64[1, 2, 3]
        b = Float64[2, 3, 5]

        expected = Point3{Float64}[[1.0, 2.0, 1.0], [2.0, 2.0, 1.0], 
        [3.0, 2.0, 1.0], [1.0, 3.0, 1.0], [2.0, 3.0, 1.0], [3.0, 3.0, 1.0], 
        [1.0, 5.0, 1.0], [2.0, 5.0, 1.0], [3.0, 5.0, 1.0], [1.0, 2.0, 2.0],
         [2.0, 2.0, 2.0], [3.0, 2.0, 2.0], [1.0, 3.0, 2.0], [2.0, 3.0, 2.0], 
         [3.0, 3.0, 2.0], [1.0, 5.0, 2.0], [2.0, 5.0, 2.0], [3.0, 5.0, 2.0], 
         [1.0, 2.0, 3.0], [2.0, 2.0, 3.0], [3.0, 2.0, 3.0], [1.0, 3.0, 3.0], 
         [2.0, 3.0, 3.0], [3.0, 3.0, 3.0], [1.0, 5.0, 3.0], [2.0, 5.0, 3.0], 
         [3.0, 5.0, 3.0]]

        result = gridpoints(a, b)

        @test result == expected
    end

    @testset "with 3 vectors" begin
        a = Float64[1, 2, 3]
        b = Float64[2, 3, 5]
        c = Float64[5, 6, 4]

        expected = Point3{Float64}[[1.0, 2.0, 5.0], [2.0, 2.0, 5.0], 
        [3.0, 2.0, 5.0], [1.0, 3.0, 5.0], [2.0, 3.0, 5.0], [3.0, 3.0, 5.0], 
        [1.0, 5.0, 5.0], [2.0, 5.0, 5.0], [3.0, 5.0, 5.0], [1.0, 2.0, 6.0], 
        [2.0, 2.0, 6.0], [3.0, 2.0, 6.0], [1.0, 3.0, 6.0], [2.0, 3.0, 6.0], 
        [3.0, 3.0, 6.0], [1.0, 5.0, 6.0], [2.0, 5.0, 6.0], [3.0, 5.0, 6.0], 
        [1.0, 2.0, 4.0], [2.0, 2.0, 4.0], [3.0, 2.0, 4.0], [1.0, 3.0, 4.0], 
        [2.0, 3.0, 4.0], [3.0, 3.0, 4.0], [1.0, 5.0, 4.0], [2.0, 5.0, 4.0], 
        [3.0, 5.0, 4.0]]

        result = gridpoints(a, b, c)

        @test result == expected
    end

    @testset "Equality of several results" begin
        a = Float64[1, 2, 3]

        result1 = gridpoints(a)
        result2 = gridpoints(a, a)
        result3 = gridpoints(a, a, a)

        @test allequal([result1, result2, result3])
    end
end


@testset "gridpoints_equilateral" verbose = true begin
    eps_level = 1e-4
       
    @testset "Vectors, non-rectangular" begin 
        xSpan = [-3,3]
        ySpan = [-2,2]
        pointSpacing = 0.5
        V = gridpoints_equilateral(xSpan,ySpan,pointSpacing; return_faces = false, rectangular=false)
        ind = round.(Int64,range(1,length(V),10))
        
        @test isa(V,Vector{Point{3,Float64}})
        @test isapprox(V[ind],Point3{Float64}[[-3.125, -2.0, 0.0], 
        [-2.375, -1.5669872981077808, 0.0], [-1.625, -1.1339745962155614, 0.0], 
        [-0.875, -0.700961894323342, 0.0], [-0.625, -0.2679491924311228, 0.0], 
        [0.625, 0.16506350946109638, 0.0], [0.875, 0.598076211353316, 0.0], 
        [1.625, 1.0310889132455352, 0.0], [2.375, 1.4641016151377544, 0.0], 
        [3.125, 1.8971143170299736, 0.0]],atol=eps_level)
        x = [v[1] for v in V]
        @test isapprox(minimum(x),xSpan[1]-pointSpacing/4,atol=eps_level)
        @test isapprox(maximum(x),xSpan[2]+pointSpacing/4,atol=eps_level)
    end

    @testset "Tuples, rectangular" begin
        xSpan = (-3,3)
        ySpan = (-2,2)
        pointSpacing = 1
        V = gridpoints_equilateral(xSpan,ySpan,pointSpacing; return_faces = false, rectangular=true)
        ind = round.(Int64,range(1,length(V),10))

        @test isa(V,Vector{Point{3,Float64}})
        @test isapprox(V[ind],Point3{Float64}[[-3.0, -2.0, 0.0], [0.75, -2.0, 0.0], 
        [-1.75, -1.1339745962155614, 0.0], [1.25, -1.1339745962155614, 0.0], 
        [-2.25, -0.2679491924311228, 0.0], [1.75, -0.2679491924311228, 0.0], 
        [-0.75, 0.598076211353316, 0.0], [2.25, 0.598076211353316, 0.0], 
        [-1.25, 1.4641016151377544, 0.0], [3.0, 1.4641016151377544, 0.0]],atol=eps_level)

        x = [v[1] for v in V]
        @test isapprox(minimum(x),xSpan[1],atol=eps_level)
        @test isapprox(maximum(x),xSpan[2],atol=eps_level)
    end


    @testset "Return faces" begin
        xSpan = (-3,3)
        ySpan = (-2,2)
        pointSpacing = 1
        F,V = gridpoints_equilateral(xSpan,ySpan,pointSpacing; return_faces = true, rectangular=true)
        ind = round.(Int64,range(1,length(V),10))
        indF = round.(Int64,range(1,length(F),10))

        @test isa(F,Vector{TriangleFace{Int64}})
        @test F[indF] == TriangleFace{Int64}[TriangleFace(1, 2, 8), 
        TriangleFace(16, 23, 22), TriangleFace(9, 10, 17), TriangleFace(3, 4, 10), 
        TriangleFace(18, 25, 24), TriangleFace(11, 12, 19), TriangleFace(33, 32, 25), 
        TriangleFace(20, 27, 26), TriangleFace(13, 14, 21), TriangleFace(35, 34, 27)]

        @test isa(V,Vector{Point{3,Float64}})
        @test isapprox(V[ind],Point3{Float64}[[-3.0, -2.0, 0.0], [0.75, -2.0, 0.0], 
        [-1.75, -1.1339745962155614, 0.0], [1.25, -1.1339745962155614, 0.0], 
        [-2.25, -0.2679491924311228, 0.0], [1.75, -0.2679491924311228, 0.0], 
        [-0.75, 0.598076211353316, 0.0], [2.25, 0.598076211353316, 0.0], 
        [-1.25, 1.4641016151377544, 0.0], [3.0, 1.4641016151377544, 0.0]],atol=eps_level)

        x = [v[1] for v in V]
        @test isapprox(minimum(x),xSpan[1],atol=eps_level)
        @test isapprox(maximum(x),xSpan[2],atol=eps_level)
    end


end


@testset "interp_biharmonic_spline" verbose = true begin
    eps_level = 1e-4

    @testset "Ranged input" begin
        x = range(0,2,12)
        y = range(0,2,12)
        xi = range(0.5,1.5,5)
        result = interp_biharmonic_spline(x, y, xi; extrapolate_method=:linear, pad_data=:linear)
        true_result = [0.5003539222091904, 0.7500704989328162, 1.0002266142134255, 
                       1.250045599106698, 1.4999641453937436]
        @test isapprox(result, true_result, atol=eps_level)
    end

    @testset "Errors" begin
        x = range(0,2,12)
        y = range(0,2,12)
        xi = range(0.5,1.5,5)
        @test_throws Exception interp_biharmonic_spline(x, y, xi; extrapolate_method=:wrong)
        @test_throws Exception interp_biharmonic_spline(x, y, xi; extrapolate_method=:linear, pad_data=:wrong)
    end

    @testset "linear interp only / linear" begin
        x = Float64[0.0, 1.0, 2.0, 3.0]
        y = Float64[0.0, 1.0, 0.0, 1.0]
        xi = range(0.5, 2.5, 5)
        result = interp_biharmonic_spline(x, y, xi; extrapolate_method=:linear)
        true_result = [0.650942317501349, 0.9999999999999994, 0.501564606542732, 
                      -3.0531133177191805e-16, 0.3537866863312682]
        @test isapprox(result, true_result, atol=eps_level)
    end

    @testset "linear interp only / linear" begin
        x = Float64[0.0, 1.0, 2.0, 3.0]
        y = Float64[0.0, 1.0, 0.0, 1.0]
        xi = range(0.5, 2.5, 5)
        result = interp_biharmonic_spline(x, y, xi; extrapolate_method=:constant)
        true_result = [0.650942317501349, 0.9999999999999994, 0.501564606542732, 
                      -3.0531133177191805e-16, 0.3537866863312682]
        @test isapprox(result, true_result, atol=eps_level)
    end

    @testset "linear / linear" begin
        x = Float64[0.0, 1.0, 2.0, 3.0]
        y = Float64[0.0, 1.0, 0.0, 1.0]
        xi = range(-0.5, 3.5, 9)
        result = interp_biharmonic_spline(x, y, xi; extrapolate_method=:linear, pad_data=:linear)
        true_result = [-0.5, -2.220446049250313e-16, 0.650942317501349,
            0.9999999999999994, 0.501564606542732,
            -2.983724378680108e-16, 0.3537866863312682,
            0.9999999999999997, 1.5]
        @test isapprox(result, true_result, atol=eps_level)
    end

    @testset "linear / constant" begin
        x = Float64[0.0, 1.0, 2.0, 3.0]
        y = Float64[0.0, 1.0, 0.0, 1.0]
        xi = range(-0.5, 3.5, 9)
        result = interp_biharmonic_spline(x, y, xi; extrapolate_method=:linear, pad_data=:constant)
        true_result = [0.0, -1.7763568394002505e-15, 0.5861167655113347,
            0.9999999999999998, 0.5015646065427324,
            -2.42861286636753e-16, 0.41861223832128147,
            0.9999999999999993, 1.0]        
        @test isapprox(result, true_result, atol=eps_level)
    end

    @testset "linear / none" begin
        x = Float64[0.0, 1.0, 2.0, 3.0]
        y = Float64[0.0, 1.0, 0.0, 1.0]
        xi = range(-0.5, 3.5, 9)
        result = interp_biharmonic_spline(x, y, xi; extrapolate_method=:linear, pad_data=:none)
        true_result = [-0.5, -1.1102230246251565e-16, 0.9548390432176067,
            0.9999999999999999, 0.5061519335211898,
            -1.1102230246251565e-16, 0.18162885699253484, 1.0, 1.5]        
        @test isapprox(result, true_result, atol=eps_level)
    end

    @testset "constant / none" begin
        x = Float64[0.0, 1.0, 2.0, 3.0]
        y = Float64[0.0, 1.0, 0.0, 1.0]
        xi = range(-0.5, 3.5, 9)
        result = interp_biharmonic_spline(x, y, xi; extrapolate_method=:constant, pad_data=:none)
        true_result = [0.0, -1.1102230246251565e-16, 0.9548390432176067,
            0.9999999999999999, 0.5061519335211898,
            -1.1102230246251565e-16, 0.18162885699253484, 1.0, 1.0]        
        @test isapprox(result, true_result, atol=eps_level)
    end

    @testset "biharmonic / none" begin
        x = Float64[0.0, 1.0, 2.0, 3.0]
        y = Float64[0.0, 1.0, 0.0, 1.0]
        xi = range(-0.5, 3.5, 9)
        result = interp_biharmonic_spline(x, y, xi; extrapolate_method=:biharmonic, pad_data=:none)
        true_result = [-2.3709643220609977, -1.1102230246251565e-16,
            0.9548390432176067, 0.9999999999999999,
            0.5061519335211898, -1.1102230246251565e-16,
            0.1816288569925348, 1.0, 2.801059658186898]        
        @test isapprox(result, true_result, atol=eps_level)
    end
end


@testset "interp_biharmonic" verbose = true begin
    eps_level = 1e-4    
    @testset "3D points 1D data, vectors" begin
        result = interp_biharmonic([[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]], [-10, 10], [[0.0, 0.0, x] for x in range(-1, 1, 5)])
        true_result = [-10.0, -7.449961786934791, 0.0, 7.449961786934791, 10.0]        
        @test isapprox(result, true_result, atol=eps_level)
    end

    @testset "3D points 1D data, geometry basics point vectors" begin
        result = interp_biharmonic(Point3{Float64}[[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]], [-10, 10],
            [Point3{Float64}(0.0, 0.0, x) for x in range(-1, 1, 5)])
        true_result = [-10.0, -7.449961786934791, 0.0, 7.449961786934791, 10.0]        
        @test isapprox(result, true_result, atol=eps_level)
    end
end


@testset "nbezier" verbose = true begin 
    eps_level = 1e-4    

    @testset "Errors" begin 
        P = Point3{Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0]]
        @test_throws Exception nbezier(P,1)
    end

    @testset "Vector Point3" begin 
        P = Point3{Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0]]
        n = 10 # Number of points
        V = nbezier(P,n) # Get Bezier fit points
        expected = Point3{Float64}[[0.0, 0.0, 0.0], 
        [0.29766803840877915, 0.034293552812071325, 0.0013717421124828531], 
        [0.5294924554183813, 0.1262002743484225, 0.010973936899862825], 
        [0.7037037037037037, 0.25925925925925924, 0.037037037037037035], 
        [0.8285322359396433, 0.4170096021947874, 0.0877914951989026], 
        [0.9122085048010974, 0.5829903978052127, 0.1714677640603567], 
        [0.9629629629629629, 0.7407407407407407, 0.2962962962962963], 
        [0.9890260631001372, 0.8737997256515775, 0.4705075445816187], 
        [0.9986282578875172, 0.9657064471879286, 0.7023319615912208], 
        [1.0, 1.0, 1.0]] 
        @test typeof(V) == typeof(P)  
        @test length(V) == n
        @test isapprox(V, expected, atol = eps_level)
    end

end 


@testset "lerp" verbose = true begin 

    @testset "Errors" begin 
        @test_throws DimensionMismatch lerp([0.0,1.0],[0.0,10.0,5.0],0.5)
    end

    @testset "1D" begin         
        @test lerp([0.0,1.0],[0.0,10.0],0.5) == 5.0 # Single value interpolation site
        @test lerp([0.0,1.0],[0.0,10.0],range(0.0,1.0,3)) == [0.0,5.0,10.0] # Range of sites
        @test lerp([0.0,1.0],[0.0,10.0],range(0.0,1.0,3)) == [0.0,5.0,10.0] # Range of sites    
        @test lerp([0.0,1.0],[0.0,10.0],[0.0,0.5,1.0]) == [0.0,5.0,10.0] # Vector of sites
        @test lerp(range(0.0,1.0,3),range(0.0,10.0,3),range(0.0,1.0,3)) == [0.0,5.0,10.0] # ranged sites, data, and values
    end

    @testset "3D points" begin 
        eps_level = 1e-4
        np = 10
        t = range(0.0,2.0*π,np) # Parameterisation metric
        V = [GeometryBasics.Point{3, Float64}(cos(t[i]),sin(t[i]),t[i]/(2.0*π)) for i in eachindex(t)] # ND data, here 3D points
        np_i = np*3 
        ti = range(minimum(t)-0.5,maximum(t)+0.5,np_i)

        @test isapprox(lerp(t,V,ti),Point3{Float64}[[1.167558325036443, -0.46036271447926463, -0.07957747154594766], [1.0833956815191297, -0.22912775185383769, -0.03960661143933059], [0.9992330380018164, 0.0021072107715893167, 0.0003642486672864905], [0.9150703944845031, 0.23334217339701627, 0.04033510877390357], [0.8309077509671899, 0.46457713602244327, 0.08030596888052065], [0.7171768097594708, 0.6710013509613107, 0.12027682898713773], [0.504069515472875, 0.7940389046839496, 0.1602476890937548], [0.29096222118627935, 0.9170764584065885, 0.2002185492003719], [0.06471611212984507, 0.9656000907937172, 0.24018940930698895], [-0.17762056150557648, 0.9228695968166506, 0.28016026941360606], [-0.419957235140998, 0.8801391028395841, 0.32013112952022316], [-0.6059298257654833, 0.7397831533655452, 0.36010198962684026], [-0.7641038558835915, 0.5512786847171849, 0.40007284973345736], [-0.9222778860016999, 0.3627742160688245, 0.4400437098400744], [-0.9396926207859084, 0.12303755372263896, 0.48001456994669145], [-0.9396926207859084, -0.12303755372263873, 0.5199854300533086], [-0.9222778860017001, -0.36277421606882426, 0.5599562901599257], [-0.7641038558835918, -0.5512786847171848, 0.5999271502665426], [-0.6059298257654835, -0.7397831533655452, 0.6398980103731597], [-0.41995723514099864, -0.880139102839584, 0.6798688704797768], [-0.17762056150557712, -0.9228695968166505, 0.7198397305863939], [0.06471611212984442, -0.9656000907937171, 0.759810590693011], [0.29096222118627946, -0.9170764584065884, 0.7997814507996283], [0.5040695154728753, -0.7940389046839496, 0.8397523109062452], [0.7171768097594708, -0.6710013509613109, 0.8797231710128623], [0.8309077509671898, -0.46457713602244327, 0.9196940311194793], [0.9150703944845031, -0.23334217339701638, 0.9596648912260964], [0.9992330380018164, -0.0021072107715894555, 0.9996357513327135], [1.0833956815191297, 0.22912775185383755, 1.0396066114393308], [1.167558325036443, 0.46036271447926436, 1.0795774715459476]],atol=eps_level)
    end

    @testset "ND" begin
        @test lerp([0.0,1.0],[[0.0,10.0,20.0,30.0],[10.0,20.0,30.0,40.0]],0.5) == [5.0,15.0,25.0,35.0] # Single value interpolation site
        @test lerp([0.0,1.0],[[0.0,10.0,20.0,30.0],[10.0,20.0,30.0,40.0]],[0.0,0.5,1.0]) == [[0.0,10.0,20.0,30.0],[5.0,15.0,25.0,35.0],[10.0,20.0,30.0,40.0]] # Single value interpolation site
    end

end

@testset "lerp_" begin 
    @test Comodo.lerp_([0.0,1.0],[0.0,10.0],0.5) == 5.0 # Vector input
    @test Comodo.lerp_(range(0.0,1.0,3),range(0.0,10.0,3),0.5) == 5.0 # range input
end


@testset "dist" verbose = true begin
    eps_level = 1e-4

    @testset "vector to vector" begin
        v1 = Float64[0, 0, 0]
        v2 = Float64[0, 0, 5]
        result = dist(v1, v2)
        @test result == [0.0 0.0 5.0; 0.0 0.0 5.0; 0.0 0.0 5.0]
    end

    @testset "vectors to vector" begin         
        v1 = [[1, 2, 3], [0, 0, 0]]
        v2 = [0, 0, 0]
        result = dist(v1, v2)
        @test result isa Matrix
        @test isapprox(result, [3.7416573867739413; 0.0;;], atol = eps_level)
    end 

    @testset "vector to vectors" begin     
        v1 = [[1, 2, 3], [0, 0, 0]]
        v2 = [0, 0, 0]
        result = dist(v2, v1)
        @test result isa Matrix
        @test isapprox(result, [3.7416573867739413 0.0], atol = eps_level)
    end 

    @testset "vector of points to vector of points" begin
        V1 = Vector{GeometryBasics.Point{3,Float64}}(undef, 4)
        V1[1] = GeometryBasics.Point{3,Float64}(1.0, 0.0, 0.0)
        V1[2] = GeometryBasics.Point{3,Float64}(0.0, 1.0, 0.0)
        V1[3] = GeometryBasics.Point{3,Float64}(0.0, 0.0, 1.0)
        V1[4] = GeometryBasics.Point{3,Float64}(1.0, 1.0, 1.0)

        V2 = Vector{GeometryBasics.Point{3,Float64}}(undef, 3)
        V2[1] = GeometryBasics.Point{3,Float64}(π, 0.0, 0.0)
        V2[2] = GeometryBasics.Point{3,Float64}(0.0, π, 0.0)
        V2[3] = GeometryBasics.Point{3,Float64}(0.0, 0.0, π)

        result = dist(V1, V2)
        eps_level = maximum(eps.(result))


        @test size(result,1) == length(V1)
        @test size(result,2) == length(V2)
        @test isapprox(result, [2.141592653589793 3.296908309475615 3.296908309475615;
                3.296908309475615 2.141592653589793 3.296908309475615;
                3.296908309475615 3.296908309475615 2.141592653589793;
                2.5664019743426345 2.5664019743426345 2.5664019743426345], atol=eps_level)
    end
end


@testset "mindist" begin     
    eps_level = 1e-4

    # Basic use
    V1 = Point3{Float64}[[1.0, 2.0, 3.0], [0.0, 0.0, 0.0], [pi, -10.0, 2.5]]
    V2 = Point3{Float64}[[4.0, 5.0, 6.0], [0.0, 0.0, 0.0]]
    D = mindist(V1, V2)    
    @test D isa Vector{Float64}
    @test isapprox(D, [3.7416573867739413, 0.0, 10.775880678677236], atol = eps_level)
    
    # Also outputting index
    D,ind = mindist(V1, V2; getIndex=true)    
    @test D isa Vector{Float64}
    @test isapprox(D, [3.7416573867739413, 0.0, 10.775880678677236], atol = eps_level)
    @test ind == [2,2,2]

    # Snap to self for self distances
    D,ind = mindist(V1, V1; getIndex=true,skipSelf = false)    
    @test D isa Vector{Float64}
    @test isapprox(D, [0.0,0.0,0.0], atol = eps_level)
    @test ind == [1,2,3]

    # Do not snap to self if self is avoided
    D,ind = mindist(V1, V1; getIndex=true,skipSelf = true)    
    @test D isa Vector{Float64}
    @test isapprox(D, [3.7416573867739413, 3.7416573867739413, 10.775880678677236], atol = eps_level)
    @test ind == [2,1,2]
end 


@testset "unique_dict_index" begin 
    result1, result2 = Comodo.unique_dict_index([1, 2, 3, 3, 3, 4, 4, 4, 5])
    @test result1 == [1, 2, 3, 4, 5]
    @test result2 == [1, 2, 3, 6, 9]

    result1, result2 = Comodo.unique_dict_index([[1, 2, 3], [3,2,1],[4,5,6]],sort_entries=true)
    @test result1 == [[1, 2, 3], [4, 5,6]]
    @test result2 == [1, 3]
end 


@testset "unique_dict_index_inverse" begin 
    result1, result2, result3 = Comodo.unique_dict_index_inverse([1, 2, 3, 3, 3, 4, 4, 4, 5])
    @test result1 == [1, 2, 3, 4, 5]
    @test result2 == [1, 2, 3, 6, 9]
    @test result3 == [1, 2, 3, 3, 3, 4, 4, 4, 5]

    result1, result2, result3 = Comodo.unique_dict_index_inverse([[1, 2, 3], [3,2,1],[4,5,6]],sort_entries=true)
    @test result1 == [[1, 2, 3], [4, 5,6]]
    @test result2 == [1, 3]
    @test result3 == [1, 1, 2]
end 


@testset "unique_dict_index_count" begin 
    result1, result2, result3 = Comodo.unique_dict_index_count([1, 2, 3, 3, 3, 4, 4, 4, 5])
    @test result1 == [1, 2, 3, 4, 5]
    @test result2 == [1, 2, 3, 6, 9]
    @test result3 == [1, 1, 3, 3, 1]

    result1, result2, result3 = Comodo.unique_dict_index_count([[1, 2, 3], [3,2,1],[4,5,6]],sort_entries=true)
    @test result1 == [[1, 2, 3], [4, 5,6]]
    @test result2 == [1, 3]
    @test result3 == [2, 1]
end 


@testset "unique_dict_index_inverse_count" begin 
    r1, r2, r3, r4 = Comodo.unique_dict_index_inverse_count([1, 2, 3, 3, 3, 4, 4, 4, 5])
    @test r1 == [1, 2, 3, 4, 5]
    @test r2 == [1, 2, 3, 6, 9]
    @test r3 == [1, 2, 3, 3, 3, 4, 4, 4, 5]
    @test r4 == [1, 1, 3, 3, 1]

    r1, r2, r3, r4 = Comodo.unique_dict_index_inverse_count([[1, 2, 3], [3,2,1],[4,5,6]],sort_entries=true)
    @test r1 == [[1, 2, 3], [4, 5,6]]
    @test r2 == [1, 3]
    @test r3 == [1, 1, 2]
    @test r4 == [2,1]
end 


@testset "unique_dict_count" begin 
    result1, result2 = Comodo.unique_dict_count([1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5])
    @test result1 == [1, 2, 3, 4, 5]
    @test result2 == [3, 4, 2, 1, 1]

    result1, result2 = Comodo.unique_dict_count([[1, 2, 3], [3,2,1],[4,5,6]],sort_entries=true)
    @test result1 == [[1, 2, 3], [4, 5,6]]
    @test result2 == [2,1]
end


@testset "unique_dict_inverse" begin 
    result1, result2 = Comodo.unique_dict_inverse([1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5])
    @test result1 == [1, 2, 3, 4, 5]
    @test result2 == [1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5]

    result1, result2 = Comodo.unique_dict_inverse([[1, 2, 3], [3,2,1],[4,5,6]],sort_entries=true)
    @test result1 == [[1, 2, 3], [4, 5,6]]
    @test result2 == [1, 1, 2]
end 

@testset "unique_dict" begin 
    result1, result2, result3 = unique_dict([1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5])
    @test result1 == [1, 2, 3, 4, 5]
    @test result2 == [1, 4, 8, 10, 11]
    @test result3 == [1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5]
end 


@testset "gunique" begin     
    r1, r2, r3, r4 = gunique([1, 2, 3, 3, 3, 4, 4, 4, 5]; return_unique=true, return_index=true, return_inverse=true, return_counts=true, sort_entries=false)
    @test r1 == [1, 2, 3, 4, 5]
    @test r2 == [1, 2, 3, 6, 9]
    @test r3 == [1, 2, 3, 3, 3, 4, 4, 4, 5]
    @test r4 == [1, 1, 3, 3, 1]

    r1, r2, r3 = gunique([1, 2, 3, 3, 3, 4, 4, 4, 5]; return_unique=true, return_index=true, return_inverse=true, return_counts=false, sort_entries=false)
    @test r1 == [1, 2, 3, 4, 5]
    @test r2 == [1, 2, 3, 6, 9]
    @test r3 == [1, 2, 3, 3, 3, 4, 4, 4, 5]

    r1, r2 = gunique([1, 2, 3, 3, 3, 4, 4, 4, 5]; return_unique=true, return_index=false, return_inverse=false, return_counts=true, sort_entries=false)
    @test r1 == [1, 2, 3, 4, 5]
    @test r2 == [1, 1, 3, 3, 1]

    r1, r2 = gunique([1, 2, 3, 3, 3, 4, 4, 4, 5]; return_unique=true, return_index=false, return_inverse=true, return_counts=false, sort_entries=false)
    @test r1 == [1, 2, 3, 4, 5]
    @test r2 == [1, 2, 3, 3, 3, 4, 4, 4, 5]

    r1, r2 = gunique([1, 2, 3, 3, 3, 4, 4, 4, 5]; return_unique=true, return_index=true, return_inverse=false, return_counts=false, sort_entries=false)
    @test r1 == [1, 2, 3, 4, 5]
    @test r2 == [1, 2, 3, 6, 9]

    r1 = gunique([1, 2, 3, 3, 3, 4, 4, 4, 5]; return_unique=true, return_index=false, return_inverse=false, return_counts=false, sort_entries=false)
    @test r1 == [1, 2, 3, 4, 5]
end 


@testset "unique_simplices" verbose = true begin
    @testset "Single triangle" begin
        F = [TriangleFace{Int64}(1, 2, 3)]       
        V = [Point3{Float64}(rand(3)) for _ in 1:3]
        F_uni, ind1, ind2 = unique_simplices(F)
        @test F_uni == F
        F_uni, ind1, ind2 = unique_simplices(F,V)
        @test F_uni == F
    end

    @testset "Set of two triangles" begin
        F = [TriangleFace{Int64}(1, 2, 3),TriangleFace{Int64}(1, 2, 3)]   
        V = [Point3{Float64}(rand(3)) for _ in 1:3]    
        F_uni, ind1, ind2 = unique_simplices(F)
        @test F_uni == [F[1]]
        @test F_uni == F[ind1]
        @test F_uni[ind2] == F
        F_uni, ind1, ind2 = unique_simplices(F,V)
        @test F_uni == [F[1]]
        @test F_uni == F[ind1]
        @test F_uni[ind2] == F
    end
end


@testset "ind2sub" verbose = true begin

    @testset "Errors" begin
        A = rand(6,8)

        # Test for out of range indices
        @test_throws BoundsError ind2sub(size(A),-1)
        @test_throws BoundsError ind2sub(size(A),length(A)+1)
    end

    ind = [1,2,3,4,8,12,30]

    @testset "1D i.e. Vector" begin
        A = rand(30)
        IJK_A = ind2sub(size(A),ind)
        @test all([A[ind[i]] == A[IJK_A[i][1]] for i in eachindex(ind)])
        @test isempty(ind2sub(size(A),Int64[])) # Check if empty is returned
    end

    @testset "2D i.e. 2D Matrix" begin
        B = rand(5,6) 
        IJK_B = ind2sub(size(B),ind)
        @test all([B[ind[i]] == B[IJK_B[i][1],IJK_B[i][2]] for i in eachindex(ind)])
    end

    @testset "3D i.e. 3D matrix" begin
        C = rand(3,5,2)
        IJK_C = ind2sub(size(C),ind)
        @test all([C[ind[i]] == C[IJK_C[i][1],IJK_C[i][2],IJK_C[i][3]] for i in eachindex(ind)])
    end

    @testset "Vector specifying size" begin
        C = rand(3,5,2)
        IJK_C = ind2sub(collect(size(C)),ind)
        @test all([C[ind[i]] == C[IJK_C[i][1],IJK_C[i][2],IJK_C[i][3]] for i in eachindex(ind)])
    end

    @testset "Tuple specifying indices" begin
        C = rand(3,5,2)
        ind_tuple = Tuple(ind[i] for i in eachindex(ind))
        IJK_C = ind2sub(collect(size(C)),ind_tuple)
        @test all([C[ind_tuple[i]] == C[IJK_C[i][1],IJK_C[i][2],IJK_C[i][3]] for i in eachindex(ind_tuple)])
    end
end

@testset "ind2sub_" verbose = true begin
    ind = 6
    @testset "1D i.e. Vector" begin
        A = rand(30)        
        @test Comodo.ind2sub_(6,length(size(A)),cumprod(size(A))) == [6]
    end

    @testset "2D i.e. 2D Matrix" begin
        B = rand(5,6)         
        @test Comodo.ind2sub_(6,length(size(B)),cumprod(size(B))) == [1,2]
    end

    @testset "3D i.e. 3D matrix" begin
        C = rand(3,5,2)        
        @test Comodo.ind2sub_(6,length(size(C)),cumprod(size(C))) == [3,2,1]
    end
end


@testset "sub2ind" verbose = true begin
    ind = [1,2,3,4,8,12,30]
    A = rand(30)
    IJK_A = [[1],[2],[3],[4],[8],[12],[30]]
    B = rand(5,6) 
    IJK_B = [[1, 1], [2, 1], [3, 1], [4, 1], [3, 2], [2, 3], [5, 6]]
    C = rand(3,5,2)
    IJK_C = [[1, 1, 1], [2, 1, 1], [3, 1, 1], [1, 2, 1], [2, 3, 1], [3, 4, 1], [3, 5, 2]]

    @testset "Errors" begin
        @test_throws DomainError sub2ind(size(A),[[-1]])
        @test_throws DomainError sub2ind(size(A),[-1]) 
        @test_throws DomainError sub2ind(size(A),[length(A)+1]) 
        @test_throws DimensionMismatch sub2ind(size(A),[1,2,3,4]) 

        @test_throws DomainError sub2ind(size(B),[[-1,1]])
        @test_throws DomainError sub2ind(size(B),[-1,1]) 
        @test_throws DomainError sub2ind(size(B),[1,length(B)+1]) 
        @test_throws DimensionMismatch sub2ind(size(B),[1,2,3,4]) 
        @test_throws DimensionMismatch sub2ind(size(B),[[1,2,3,4]])

        @test_throws DomainError sub2ind(size(C),[[-1,1,1]])
        @test_throws DomainError sub2ind(size(C),[-1,1,1]) 
        @test_throws DomainError sub2ind(size(C),[length(C)+1,1,1]) 
        @test_throws DimensionMismatch sub2ind(size(C),[1,2,3,4]) 
        @test_throws DimensionMismatch sub2ind(size(B),[[1,2,3,4]])
    end

    @testset "1D i.e. Vector" begin        
        @test sub2ind(size(A),IJK_A)==ind
    end

    @testset "2D i.e. 2D Matrix" begin    
        @test sub2ind(size(B),IJK_B)==ind
    end

    @testset "3D i.e. 3D matrix" begin        
        @test sub2ind(size(C),IJK_C)==ind
    end

    @testset "Vector specifying indices 1D" begin        
        @test [sub2ind(size(A),[i])[1] for i in ind]==ind # IJK = ind for 1D case
    end

    @testset "Vector specifying size" begin        
        @test sub2ind(collect(size(C)),IJK_C)==ind
    end
end


@testset "sub2ind_" verbose = true begin

    @testset "1D i.e. Vector" begin
        A = rand(30)        
        @test Comodo.sub2ind_([6],length(size(A)),cumprod(size(A))) == 6
    end

    @testset "2D i.e. 2D Matrix" begin
        B = rand(5,6)         
        @test Comodo.sub2ind_([1,2],length(size(B)),cumprod(size(B))) == 6
    end

    @testset "3D i.e. 3D matrix" begin
        C = rand(3,5,2)        
        @test Comodo.sub2ind_([3,2,1],length(size(C)),cumprod(size(C))) == 6
    end
end


@testset "meshedges" verbose = true begin
    @testset "Single triangle" begin
        F = [TriangleFace{Int64}(1, 2, 3)]       
        E = meshedges(F)
        @test E == LineFace{Int64}[[1, 2], [2, 3], [3, 1]]
    end

    @testset "Single quad" begin
        F = [QuadFace{Int64}(1, 2, 3, 4)]       
        E = meshedges(F)
        @test E == LineFace{Int64}[[1, 2], [2, 3], [3, 4], [4, 1]]
    end

    @testset "Triangles" begin
        F = [TriangleFace{Int64}(1, 2, 3),TriangleFace{Int64}(1, 4, 3)]
        E = meshedges(F)
        @test E == LineFace{Int64}[[1, 2], [1, 4], [2, 3], [4, 3], [3, 1], [3, 1]]
    end

    @testset "Quads" begin
        F = [QuadFace{Int64}(1, 2, 3, 4),QuadFace{Int64}(6, 5, 4, 3)]
        E = meshedges(F)
        @test E == LineFace{Int64}[[1, 2], [6, 5], [2, 3], [5, 4], [3, 4], [4, 3], [4, 1], [3, 6]]
    end

    @testset "Mesh" begin
        M = cube(1.0)
        E = meshedges(M,unique_only=true)
        @test E == LineFace{Int64}[[1, 2], [7, 8], [5, 6], [6, 7], [5, 8], [2, 3], 
        [2, 6], [3, 7], [4, 8], [1, 5], [3, 4], [1, 4]]
    end
end


@testset "icosahedron" begin
    eps_level = 1e-4
    r = 1.0
    ϕ = Base.MathConstants.golden # (1.0+sqrt(5.0))/2.0, Golden ratio
    s = r/sqrt(ϕ + 2.0)
    t = ϕ*s
    M = icosahedron(r)
    F = faces(M)
    V = coordinates(M)
    @test M isa GeometryBasics.Mesh{3,Float64,GeometryBasics.Ngon{3,Float64,3,Point3{Float64}},SimpleFaceView{3,Float64,3,Int64,Point3{Float64},TriangleFace{Int64}}}
    @test length(F) == 20
    @test isapprox(V[1], [0.0, -s, -t], atol=eps_level)
end


@testset "octahedron" begin
    eps_level = 1e-4
    r = 1.0
    s = r/sqrt(2.0)
    M = octahedron(1.0) 
    F = faces(M)
    V = coordinates(M)
    @test M isa GeometryBasics.Mesh{3,Float64,GeometryBasics.Ngon{3,Float64,3,Point3{Float64}},SimpleFaceView{3,Float64,3,Int64,Point3{Float64},TriangleFace{Int64}}}
    @test length(F) == 8
    @test isapprox(V[1], [-s,  -s, 0.0], atol=eps_level)
end


@testset "dodecahedron" begin
    eps_level = 1e-4
    r = 1.0
    ϕ = Base.MathConstants.golden # (1.0+sqrt(5.0))/2.0, Golden ratio
    s = r/sqrt(3.0)
    t = ϕ*s    
    w = (ϕ-1.0)*s
    M = dodecahedron(r)
    F = faces(M)
    V = coordinates(M)
    @test M isa GeometryBasics.Mesh{3, Float64, GeometryBasics.Ngon{3, Float64, 5, Point3{Float64}}, SimpleFaceView{3, Float64, 5, Int64, Point3{Float64}, NgonFace{5, Int64}}}
    @test length(F) == 12
    @test isapprox(V[1], [s,s,s], atol=eps_level)
end


@testset "cube" begin
    eps_level = 1e-4
    r = 1.0
    s = r/sqrt(3.0)
    M = cube(1.0) 
    F = faces(M)
    V = coordinates(M)
    @test M isa GeometryBasics.Mesh{3, Float64, GeometryBasics.Ngon{3, Float64, 4, Point3{Float64}}, SimpleFaceView{3, Float64, 4, Int64, Point3{Float64}, QuadFace{Int64}}}
    @test length(F) == 6
    @test isapprox(V[1], [-s,  -s, -s], atol=eps_level)
end


@testset "tetrahedron" begin
    eps_level = 1e-4
    r = 1.0
    a = r*sqrt(2.0)/sqrt(3.0)
    b = -r*sqrt(2.0)/3.0
    c = -r/3.0    
    M = tetrahedron(1.0) 
    F = faces(M)
    V = coordinates(M)
    @test M isa GeometryBasics.Mesh{3,Float64,GeometryBasics.Ngon{3,Float64,3,Point3{Float64}},SimpleFaceView{3,Float64,3,Int64,Point3{Float64},TriangleFace{Int64}}}
    @test length(F) == 4
    @test isapprox(V[1], [-a,  b, c], atol=eps_level)
end


@testset "platonicsolid" begin
    eps_level = 1e-4

    r = 2.5 # Radius
    lv = [4,8,6,12,20] # Correct vertex numbers
    lf = [4,6,8,20,12] # Correct vertex numbers

    for q=1:5
        M = platonicsolid(q, r) # icosahedron
        F = faces(M)
        V = coordinates(M)
        @test isa(M,GeometryBasics.Mesh)
        @test length(F) == lf[q]
        @test length(V) == lv[q]
        @test isapprox(mean(norm.(V)), r, atol=eps_level)
    end
end


@testset "tofaces" verbose = true begin
    # Edges matrix and vector
    Fem = [1 2; 4 5]
    Fev = [[1,2],[4,5]]

    # Triangles matrix and vector
    Ftm = [1 2 3; 4 5 6]
    Ftv = [[1,2,3],[4,5,6]]
    
    # Quads matrix and vector
    Fqm = [1 2 3 4; 5 6 7 8]
    Fqv = [[1,2,3,4],[5,6,7,8]]

    # Ngons (pentagons) matrix and vector
    Fnm = [1 2 3 4 5; 6 7 8 9 10]
    Fnv = [[1,2,3,4,5],[6,7,8,9,10]]

    # Imported triangular mesh 
    fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
    M = load(fileName_mesh)   
    
    @testset "Matrix input" verbose = true begin        
        Fem_G = tofaces([Fem[1,:]])
        @testset "1 LineFace" begin        
            @test isa(Fem_G,Vector{GeometryBasics.LineFace{Int64}})
            @test length(Fem_G) == 1
        end

        Ftm_G = tofaces([Ftm[1,:]])
        @testset "1 TriangleFace" begin        
            @test isa(Ftm_G,Vector{GeometryBasics.TriangleFace{Int64}})
            @test length(Ftm_G) == 1
        end

        Fqm_G = tofaces([Fqm[1,:]])
        @testset "1 QuadFace" begin        
            @test isa(Fqm_G,Vector{GeometryBasics.QuadFace{Int64}})
            @test length(Fqm_G) == 1
        end

        Fnm_G = tofaces([Fnm[1,:]])
        @testset "1 NgonFace" begin        
            @test isa(Fnm_G,Vector{GeometryBasics.NgonFace{5,Int64}})
            @test length(Fnm_G) == 1
        end  
        
        Fem_G = tofaces(Fem)
        @testset "LineFace" begin        
            @test isa(Fem_G,Vector{GeometryBasics.LineFace{Int64}})
            @test length(Fem_G) == size(Fem,1)
        end

        Ftm_G = tofaces(Ftm)
        @testset "TriangleFace" begin        
            @test isa(Ftm_G,Vector{GeometryBasics.TriangleFace{Int64}})
            @test length(Ftm_G) == size(Ftm,1)
        end

        Fqm_G = tofaces(Fqm)
        @testset "QuadFace" begin        
            @test isa(Fqm_G,Vector{GeometryBasics.QuadFace{Int64}})
            @test length(Fqm_G) == size(Fqm,1)
        end

        Fnm_G = tofaces(Fnm)
        @testset "NgonFace" begin        
            @test isa(Fnm_G,Vector{GeometryBasics.NgonFace{5,Int64}})
            @test length(Fnm_G) == size(Fnm,1)
        end        
    end

    @testset "vector input" verbose = true begin        
        Fev_G = tofaces([Fev[1]])        
        @testset "1 LineFace" begin        
            @test isa(Fev_G,Vector{GeometryBasics.LineFace{Int64}})
            @test length(Fev_G) == 1
        end

        Ftv_G = tofaces([Ftv[1]])        
        @testset "1 TriangleFace" begin        
            @test isa(Ftv_G,Vector{GeometryBasics.TriangleFace{Int64}})
            @test length(Ftv_G) == 1
        end

        Fqv_G = tofaces([Fqv[1]])
        @testset "1 QuadFace" begin        
            @test isa(Fqv_G,Vector{GeometryBasics.QuadFace{Int64}})
            @test length(Fqv_G) == 1
        end

        Fnv_G = tofaces([Fnv[1]])
        @testset "1 NgonFace" begin        
            @test isa(Fnv_G,Vector{GeometryBasics.NgonFace{5,Int64}})
            @test length(Fnv_G) == 1
        end  

        Fev_G = tofaces(Fev)        
        @testset "LineFace" begin        
            @test isa(Fev_G,Vector{GeometryBasics.LineFace{Int64}})
            @test length(Fev_G) == length(Fev)
        end

        Ftv_G = tofaces(Ftv)
        @testset "TriangleFace" begin        
            @test isa(Ftv_G,Vector{GeometryBasics.TriangleFace{Int64}})
            @test length(Ftv_G) == length(Ftv)
        end

        Fqv_G = tofaces(Fqv)
        @testset "QuadFace" begin        
            @test isa(Fqv_G,Vector{GeometryBasics.QuadFace{Int64}})
            @test length(Fqv_G) == length(Fqv)
        end

        Fnv_G = tofaces(Fnv)
        @testset "NgonFace" begin        
            @test isa(Fnv_G,Vector{GeometryBasics.NgonFace{5,Int64}})
            @test length(Fnv_G) == length(Fnv)
        end    
        
        Fnv_G2 = tofaces(Fnv_G)
        @testset "NgonFace no change" begin        
            @test isa(Fnv_G2,Vector{GeometryBasics.NgonFace{5,Int64}})
            @test length(Fnv_G2) == length(Fnv_G)
        end    

        @testset "Imported mesh points offset integer based" verbose = true begin
            Fme = [NgonFace{2, OffsetInteger{-1, UInt32}}(1,2),
                   NgonFace{2, OffsetInteger{-1, UInt32}}(2,3)] 
           
            Fmt = [NgonFace{3, OffsetInteger{-1, UInt32}}(1,2,3),
                   NgonFace{3, OffsetInteger{-1, UInt32}}(2,3,4)] 
            
            Fmq = [NgonFace{4, OffsetInteger{-1, UInt32}}(1,2,3,4),
                   NgonFace{4, OffsetInteger{-1, UInt32}}(3,4,5,6)]
            
            Fmn = [NgonFace{5, OffsetInteger{-1, UInt32}}(1,2,3,4,5),
                   NgonFace{5, OffsetInteger{-1, UInt32}}(3,4,5,6,7)]

            @testset "edges" begin
                F = tofaces(Fme)
                @test isa(F,Vector{GeometryBasics.LineFace{Int64}})
                @test length(F) == length(Fme)
                @test Fme[1] == F[1]
            end  

            @testset "Triangles" begin
                F = tofaces(Fmt)
                @test isa(F,Vector{GeometryBasics.TriangleFace{Int64}})
                @test length(F) == length(Fmt)
                @test Fmt[1] == F[1]
            end

            @testset "Quads" begin
                F = tofaces(Fmq)
                @test isa(F,Vector{GeometryBasics.QuadFace{Int64}})
                @test length(F) == length(Fmq)
                @test Fmq[1] == F[1]
            end

            @testset "NgonFace" begin
                F = tofaces(Fmn)
                @test isa(F,Vector{GeometryBasics.NgonFace{5,Int64}})
                @test length(F) == length(Fmn)
                @test Fmn[1] == F[1]
            end
        end
    end
end


@testset "topoints" verbose = true begin

    @testset "Matrix input" begin        
        V = topoints(rand(10,3))
        @test isa(V,Vector{GeometryBasics.Point3{Float64}})
        @test length(V) == 10
    end

    @testset "Vector Float64" begin
        Vv = [rand(3) for _ in 1:5]       
        V = topoints(Vv)
        @test isa(V,Vector{GeometryBasics.Point3{Float64}})
        @test length(V) == 5
    end

    @testset "Vector Vec3" begin
        Vv = Vector{Vec3{Float64}}(undef,5)       
        V = topoints(Vv)
        @test isa(V,Vector{GeometryBasics.Point3{Float64}})
        @test length(V) == 5
    end

    @testset "Vector Vec{m,Float64}}" begin
        m = 4
        Vv = Vector{Vec{m,Float64}}(undef,5)       
        V = topoints(Vv)
        @test isa(V,Vector{GeometryBasics.Point{m,Float64}})
        @test length(V) == 5

        m = 5
        Vv = Vector{Vec{m,Float64}}(undef,5)       
        V = topoints(Vv)
        @test isa(V,Vector{GeometryBasics.Point{m,Float64}})
        @test length(V) == 5
    end

    @testset "Imported mesh points" begin
        # Imported triangular mesh 
        fileName_mesh = joinpath(comododir(),"assets","stl","stanford_bunny_low.stl")
        M = load(fileName_mesh) 
        Vv = coordinates(M)       
        V = topoints(Vv)
        @test isa(V,Vector{GeometryBasics.Point3{Float32}})
        @test length(V) == length(Vv)
    end

    @testset "points no change" begin
        # Imported triangular mesh 
        Vv = rand(Point{3,Float64},5) 
        V = topoints(Vv)
        @test isa(V,Vector{GeometryBasics.Point3{Float64}})
        @test length(V) == length(Vv)
    end

end


@testset "togeometrybasics_mesh" verbose = true begin

    @testset "Matrix input" begin  
        @testset "triangles" begin      
            n = 10 
            m = 3
            k = 5
            M = togeometrybasics_mesh(rand(n,3),rand(1:n,k,m))
            @test isa(M,GeometryBasics.Mesh)
            @test length(coordinates(M))==n
            @test length(faces(M))==k
            @test length(faces(M)[1])==m
        end
        @testset "quads" begin      
            n = 12 
            m = 4
            k = 6
            M = togeometrybasics_mesh(rand(n,3),rand(1:n,k,m))
            @test isa(M,GeometryBasics.Mesh)
            @test length(coordinates(M))==n
            @test length(faces(M))==k
            @test length(faces(M)[1])==m
        end
    end

    @testset "Vector input" begin  
        @testset "triangles" begin      
            n = 10 
            m = 3
            k = 5
            V = Vector{Vec3{Float64}}(undef,n)
            F = [rand(1:n,m) for i=1:k]
            M = togeometrybasics_mesh(V,F)
            @test isa(M,GeometryBasics.Mesh)
            @test length(coordinates(M))==n
            @test length(faces(M))==k
            @test length(faces(M)[1])==m
        end
        @testset "quads" begin      
            n = 12 
            m = 4
            k = 6
            V = Vector{Vec3{Float64}}(undef,n)
            F = [rand(1:n,m) for i=1:k]
            M = togeometrybasics_mesh(V,F)
            @test isa(M,GeometryBasics.Mesh)
            @test length(coordinates(M))==n
            @test length(faces(M))==k
            @test length(faces(M)[1])==m
        end
    end
end

@testset "edgecrossproduct" verbose = true begin
    @testset "Single triangle" begin
        F = [TriangleFace{Int64}(1, 2, 3)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,3)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        C = edgecrossproduct(F,V) 
        @test C == [Vec3{Float64}(0.0,0.0,0.5)]
    end
    @testset "Single quad" begin
        F = [QuadFace{Int64}(1, 2, 3, 4)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        C = edgecrossproduct(F,V) 
        @test C == [Vec3{Float64}(0.0,0.0,1.0)]
    end
    @testset "Triangles" begin
        F = [TriangleFace{Int64}(1, 2, 3),TriangleFace{Int64}(1, 4, 3)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        C = edgecrossproduct(F,V) 
        @test C == [Vec3{Float64}(0.0,0.0,0.5),Vec3{Float64}(0.0,0.0,-0.5)]
    end
    @testset "Quads" begin
        F = [QuadFace{Int64}(1, 2, 3, 4),QuadFace{Int64}(6, 5, 4, 3)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,6)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        V[5] = GeometryBasics.Point{3, Float64}(2.0, 0.0, 0.0)
        V[6] = GeometryBasics.Point{3, Float64}(2.0, 1.0, 0.0)
        C = edgecrossproduct(F,V) 
        @test C == [Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,-1.0)]
    end
    @testset "Mesh" begin
        M = cube(1.0)
        C = edgecrossproduct(M) 
        @test C == Vec3{Float64}[[0.0, 0.0, -1.3333333333333337], 
        [0.0, 0.0, 1.3333333333333337], [-1.3333333333333337, 0.0, 0.0], 
        [0.0, 1.3333333333333337, 0.0], [1.3333333333333337, 0.0, 0.0], 
        [0.0, -1.3333333333333337, 0.0]]
    end
end

@testset "facenormal" verbose = true begin
    @testset "Single triangle" begin
        F = [TriangleFace{Int64}(1, 2, 3)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,3)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        N = facenormal(F,V) 
        @test N == [Vec3{Float64}(0.0,0.0,1.0)]
    end
    @testset "Single quad" begin
        F = [QuadFace{Int64}(1, 2, 3, 4)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        N = facenormal(F,V) 
        @test N == [Vec3{Float64}(0.0,0.0,1.0)]
    end
    @testset "Triangles" begin
        F = [TriangleFace{Int64}(1, 2, 3),TriangleFace{Int64}(1, 4, 3)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        N = facenormal(F,V) 
        @test N == [Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,-1.0)]
    end
    @testset "Quads" begin
        F = [QuadFace{Int64}(1, 2, 3, 4),QuadFace{Int64}(6, 5, 4, 3)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,6)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        V[5] = GeometryBasics.Point{3, Float64}(2.0, 0.0, 0.0)
        V[6] = GeometryBasics.Point{3, Float64}(2.0, 1.0, 0.0)
        N = facenormal(F,V) 
        @test N == [Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,-1.0)]
    end
    @testset "Mesh" begin
        M = cube(1)
        @test facenormal(M) == Vec3{Float64}[[0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]
    end
end

@testset "facearea" verbose = true begin
    @testset "Single triangle" begin
        F = [TriangleFace{Int64}(1, 2, 3)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,3)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        A = facearea(F,V) 
        @test A == [0.5]
    end
    @testset "Single quad" begin
        F = [QuadFace{Int64}(1, 2, 3, 4)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        A = facearea(F,V) 
        @test A == [1.0]
    end
    @testset "Triangles" begin
        F = [TriangleFace{Int64}(1, 2, 3),TriangleFace{Int64}(1, 4, 3)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        A = facearea(F,V) 
        @test A == [0.5,0.5]
    end
    @testset "Quads" begin
        F = [QuadFace{Int64}(1, 2, 3, 4),QuadFace{Int64}(6, 5, 4, 3)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,6)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        V[5] = GeometryBasics.Point{3, Float64}(2.0, 0.0, 0.0)
        V[6] = GeometryBasics.Point{3, Float64}(2.0, 1.0, 0.0)
        A = facearea(F,V) 
        @test A == [1.0,1.0]
    end
    @testset "Mesh" begin
        M = cube(1)
        @test facearea(M) == [1.3333333333333337, 1.3333333333333337, 1.3333333333333337, 1.3333333333333337, 1.3333333333333337, 1.3333333333333337]
    end
end

@testset "vertexnormal" verbose = true begin
    @testset "Single triangle" begin
        F = [TriangleFace{Int64}(1, 2, 3)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,3)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        N = vertexnormal(F,V) 
        @test N == [Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0)]
    end
    @testset "Single quad" begin
        F = [QuadFace{Int64}(1, 2, 3, 4)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        N = vertexnormal(F,V) 
        @test N == [Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0)]
    end
    @testset "Triangles" begin
        F = [TriangleFace{Int64}(1, 2, 3),TriangleFace{Int64}(1, 4, 3)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        N = vertexnormal(F,V) 
        @test N == [Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0)]
    end
    @testset "Quads" begin
        F = [QuadFace{Int64}(1, 2, 3, 4),QuadFace{Int64}(3,4,5,6)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,6)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        V[5] = GeometryBasics.Point{3, Float64}(2.0, 0.0, 0.0)
        V[6] = GeometryBasics.Point{3, Float64}(2.0, 1.0, 0.0)
        N = vertexnormal(F,V) 
        @test N == [Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0),Vec3{Float64}(0.0,0.0,1.0)]
    end
    @testset "Mesh" begin
        M = cube(1)
        @test vertexnormal(M) == Vec3{Float64}[[-0.5773502691896257, -0.5773502691896257, -0.5773502691896257], [-0.5773502691896257, 0.5773502691896257, -0.5773502691896257], [0.5773502691896257, 0.5773502691896257, -0.5773502691896257], [0.5773502691896257, -0.5773502691896257, -0.5773502691896257], [-0.5773502691896257, -0.5773502691896257, 0.5773502691896257], [-0.5773502691896257, 0.5773502691896257, 0.5773502691896257], [0.5773502691896257, 0.5773502691896257, 0.5773502691896257], [0.5773502691896257, -0.5773502691896257, 0.5773502691896257]]
    end
end

@testset "edgelengths" begin
    @testset "GeometryBasics faces, vertices" begin
        F = [QuadFace{Int64}(1, 2, 3, 4)]
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        
        @test edgelengths(F,V) == [1.0, 1.0, 1.0, 1.0] # Unit square
        @test edgelengths(F,V*pi) == pi.*[1.0, 1.0, 1.0, 1.0] # Scaled square
    end

    @testset "GeometryBasics LineFace edges" begin
        F = [QuadFace{Int64}(1, 2, 3, 4)]
        E = meshedges(F)
        V = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
        V[1] = GeometryBasics.Point{3, Float64}(0.0, 0.0, 0.0)
        V[2] = GeometryBasics.Point{3, Float64}(1.0, 0.0, 0.0)
        V[3] = GeometryBasics.Point{3, Float64}(1.0, 1.0, 0.0)
        V[4] = GeometryBasics.Point{3, Float64}(0.0, 1.0, 0.0)
        
        @test edgelengths(E,V) == [1.0, 1.0, 1.0, 1.0]        
    end

    @testset "Mesh" begin
        M = cube(sqrt(3))
        @test edgelengths(M) == 2.0*ones(12)
    end
end


@testset "subtri" verbose = true begin
    eps_level = 1e-4
    M = platonicsolid(4, 1.0) # icosahedron with radius 1.0 
    V = coordinates(M)
    F = faces(M)
    n = 2

    @testset "Errors" begin
        @test_throws ArgumentError subtri(F, V, n; method=:wrong)        
        @test_throws ArgumentError subtri(F, V, -1)
    end

    @testset "linear" begin
        Fn, Vn = subtri(F, V, 0)
        @test Fn == F
        @test Vn == V

        Fn, Vn = subtri(F, V, n; method=:linear)
        ind = round.(Int64,range(1,length(Vn),6)) # Sample indices

        @test eltype(F) == eltype(Fn)
        @test length(Fn) == 4^n*length(F)        
        @test eltype(V) == eltype(Vn)        
        @test isapprox(Vn[ind], Point{3, Float64}[[0.0, -0.5257311121191336, -0.85065080835204], 
        [-0.2628655560595668, 0.6881909602355868, 0.42532540417602], [0.0, -0.6881909602355868, -0.42532540417602], 
        [0.1314327780297834, -0.7694208842938134, 0.21266270208801], [-0.21266270208801, 0.3942983340893502, 0.7694208842938134], 
        [-0.63798810626403, 0.1314327780297834, -0.6069610361773602]], atol=eps_level)        
    end

    @testset "Loop" begin
        Fn, Vn = subtri(F, V, 0; method=:Loop)
        @test Fn == F
        @test Vn == V

        Fn, Vn = subtri(F, V, n; method=:Loop)
        ind = round.(Int64,range(1,length(Vn),6)) # Sample indices

        @test eltype(F) == eltype(Fn)
        @test length(Fn) == 4^n*length(F)        
        @test eltype(V) == eltype(Vn)        
        @test isapprox(Vn[ind], Point{3, Float64}[[-4.668111393312408e-18, -0.37854357368761526, -0.6124963684494119], 
        [-0.22191241796362113, 0.5809742527544328, 0.3590618347908117], [0.0, -0.6134756030005282, -0.37014980570299627], 
        [0.10988309631672716, -0.681387091313344, 0.19235522107345335], [-0.17398693193931355, 0.3182970619500977, 0.6225453022912995], 
        [-0.5150154647544891, 0.10752983753681045, -0.4922839938894113]], atol=eps_level)
    end

    @testset "Loop, constrained boundary" begin
        # Example with boundary edges (extruded prism)
        r = 1.0
        nc = 3
        Vc = circlepoints(r,nc;dir=:cw)    
        d = norm(Vc[1]-Vc[2])        
        s = 1    
        F,V = extrudecurve(Vc,d;s=s, num_steps=2, close_loop=true,face_type=:tri_slash)
        
        Fn, Vn = subtri(F, V, 0; method=:Loop, constrain_boundary=true)
        @test Fn == F
        @test Vn == V

        Fn, Vn = subtri(F, V, n; method=:Loop, constrain_boundary=true)
        ind = round.(Int64,range(1,length(Vn),6)) # Sample indices

        @test eltype(F) == eltype(Fn)
        @test length(Fn) == 4^n*length(F)        
        @test eltype(V) == eltype(Vn)        
        @test isapprox(Vn[ind], Point{3, Float64}[[1.0, 0.0, 0.0], [-0.4999999999999999, 1.1102230246251565e-16, 1.7320508075688776], 
        [-0.054687500000000125, -0.5277342304311423, 0.8660254037844388], [0.2734375000000002, 0.39241776108982385, 1.2990381056766582], 
        [-0.3359374999999997, 0.5818608181676699, 1.2990381056766582], [-0.1249999999999997, 0.6495190528383291, 1.7320508075688776]], atol=eps_level)
    end

end


@testset "subquad" verbose = true begin
    eps_level = 1e-4
    M = cube(1.0)
    F = faces(M)
    V = coordinates(M)    
    n = 3

    @testset "Errors" begin
        @test_throws ArgumentError subquad(F, V, n; method=:wrong)
        @test_throws ArgumentError subquad(F, V, -1)
    end

    @testset "linear" begin
        Fn, Vn = subquad(F, V, 0; method=:Linear)
        @test Fn == F
        @test Vn == V

        Fn, Vn = subquad(F, V, n; method=:linear)        
        ind = round.(Int64,range(1,length(Vn),6)) # Sample indices
        @test eltype(F) == eltype(Fn)
        @test length(Fn) == 4^n*length(F)        
        @test eltype(V) == eltype(Vn)  
        @test isapprox(Vn[ind], Point{3, Float64}[[-0.5773502691896258, -0.5773502691896258, -0.5773502691896258], 
        [0.2886751345948129, -0.2886751345948129, -0.5773502691896258], [0.4330127018922194, 0.5773502691896258, -0.5773502691896258], 
        [-0.5773502691896258, 0.2886751345948129, 0.14433756729740646], [0.14433756729740646, -0.14433756729740646, 0.5773502691896258], 
        [0.14433756729740646, -0.5773502691896258, -0.43301270189221935]], atol=eps_level)        
    end

    @testset "Catmull_Clark" begin
        Fn, Vn = subquad(F, V, 0; method=:Catmull_Clark)
        @test Fn == F
        @test Vn == V

        Fn, Vn = subquad(F, V, n; method=:Catmull_Clark)
        ind = round.(Int64,range(1,length(Vn),6)) # Sample indices
        @test eltype(F) == eltype(Fn)
        @test length(Fn) == 4^n*length(F)        
        @test eltype(V) == eltype(Vn)  
        @test isapprox(Vn[ind], Point{3, Float64}[[-0.2895661072324513, -0.2895661072324513, -0.2895661072324513], 
        [0.18414681640860342, -0.18414681640860342, -0.4257509268592806], [0.24244200756986245, 0.31179169810728824, -0.31179169810728824], 
        [-0.44889359308574917, 0.19141642225574024, 0.09422035643025145], [0.0977285611909523, -0.0977285611909523, 0.47360764269461497], 
        [0.0907121516695506, -0.407828803431474, -0.2770228830681994]], atol=eps_level)        
    end

    @testset "Catmull_Clark, constrained boundary" begin
        # Example with boundary edges (extruded prism)
        r = 1.0
        nc = 3
        Vc = circlepoints(r,nc;dir=:cw)    
        d = norm(Vc[1]-Vc[2])        
        s = 1    
        F,V = extrudecurve(Vc,d;s=s, num_steps=2, close_loop=true,face_type=:quad)

        Fn, Vn = subquad(F, V, 0; method=:Catmull_Clark, constrain_boundary=true)
        @test Fn == F
        @test Vn == V

        Fn, Vn = subquad(F, V, n; method=:Catmull_Clark, constrain_boundary=true)
        ind = round.(Int64,range(1,length(Vn),6)) # Sample indices
        @test eltype(F) == eltype(Fn)
        @test length(Fn) == 4^n*length(F)        
        @test eltype(V) == eltype(Vn)  
        @test isapprox(Vn[ind], Point{3, Float64}[[1.0, 0.0, 0.0], [0.625, -0.21650635094610965, 0.0], 
        [-0.4516601562499999, 0.16068830734281586, 0.8660254037844388], [-0.4296874999999998, 0.31757083898540706, 1.0825317547305486], 
        [-0.25976562500000017, -0.6664336115059938, 1.515544456622768], [0.07128906250000025, 0.515894039363777, 0.2165063509461097]], atol=eps_level)        
    end
end


@testset "geosphere" begin    
    eps_level = 1e-4
    r = 1.0
    
    n = 3
    F, V = geosphere(n, r)    
    ind = round.(Int64,range(1,length(V),6)) # Sample indices
    @test F isa Vector{TriangleFace{Int64}}
    @test length(F) == 4^n*20    
    @test V isa Vector{Point3{Float64}}    
    @test isapprox(V[ind],Point{3, Float64}[[0.0, -0.5257311121191336, -0.85065080835204], 
    [-0.5877852522924731, -0.6881909602355868, -0.42532540417602], [0.26286555605956685, -0.16245984811645317, -0.9510565162951538], 
    [-0.25989191300775444, 0.43388856455269487, 0.8626684804161864], [0.9129824929322992, 0.3996070517018512, 0.08232358003196016], 
    [0.13279247682790243, -0.2201170274732924, 0.9663925974024391]], atol=eps_level)
    @test isapprox(norm.(V),fill(r,length(V))) # Test radii
end


@testset "hexbox" verbose = true begin
    @testset "Single hex box" begin
        E,V,F,Fb,CFb_type = hexbox([1.0,1.0,1.0],[1,1,1])
        @test E == [[1, 2, 4, 3, 5, 6, 8, 7]] 
        @test V == Point3{Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
        @test F == QuadFace{Int64}[[3, 4, 2, 1], [5, 6, 8, 7], [1, 2, 6, 5], [4, 3, 7, 8], [2, 4, 8, 6], [3, 1, 5, 7]]
        @test Fb == QuadFace{Int64}[[3, 4, 2, 1], [5, 6, 8, 7], [1, 2, 6, 5], [4, 3, 7, 8], [2, 4, 8, 6], [3, 1, 5, 7]]
        @test CFb_type == [1, 2, 3, 4, 5, 6]
    end
    @testset "2x2x2 hex box" begin
        E,V,F,Fb,CFb_type = hexbox([1.0,1.0,1.0],[2,2,2])
        @test E[1] == [1, 2, 5, 4, 10, 11, 14, 13]
        @test length(E) == 8
        @test V[5] == [0.5, 0.5, 0.0]
        @test length(V) == 27
        @test length(F) == 48
        @test length(Fb) == 24
        @test CFb_type == [1, 3, 6, 1, 3, 5, 1, 4, 6, 1, 4, 5, 2, 3, 6, 2, 3, 5, 2, 4, 6, 2, 4, 5]
    end
end


@testset "con_face_edge" verbose = true begin    
    @testset "Single triangle" begin
        F = TriangleFace{Int64}[[1,2,3]]
        E = meshedges(F;unique_only=false)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    

        con_F2E = con_face_edge(F)
        @test con_F2E == [[1,2,3]]

        con_F2E_2 = con_face_edge(F,E_uni,indReverse)
        @test con_F2E == con_F2E_2
    end

    @testset "Single quad" begin
        F = QuadFace{Int64}[[1,2,3,4]]
        E = meshedges(F;unique_only=false)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    

        con_F2E = con_face_edge(F)
        @test con_F2E == [[1,2,3,4]]

        con_F2E_2 = con_face_edge(F,E_uni,indReverse)
        @test con_F2E == con_F2E_2
    end

    @testset "Triangles" begin
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]        
        con_F2E = con_face_edge(F)
        @test con_F2E == [[1,2,4],[2,3,5]]        
    end

    @testset "Quads" begin
        F = QuadFace{Int64}[[1,2,3,4],[3,4,5,6]]        
        con_F2E = con_face_edge(F)
        @test con_F2E == [[1,3,2,6],[2,4,5,7]]        
    end
end


@testset "con_edge_face" verbose = true begin    
    @testset "Single triangle" begin
        F = TriangleFace{Int64}[[1,2,3]]
        E = meshedges(F;unique_only=false)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    

        con_E2F = con_edge_face(F)
        @test con_E2F == fill([1],length(F[1]))

        con_E2F_2 = con_edge_face(F,E_uni,indReverse)
        @test con_E2F == con_E2F_2
    end

    @testset "Single quad" begin
        F = QuadFace{Int64}[[1,2,3,4]]
        E = meshedges(F;unique_only=false)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    

        con_E2F = con_edge_face(F)
        @test con_E2F == fill([1],length(F[1]))

        con_E2F_2 = con_edge_face(F,E_uni,indReverse)
        @test con_E2F == con_E2F_2
    end

    @testset "Triangles" begin
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]        
        con_E2F = con_edge_face(F)
        @test con_E2F == [[1], [1, 2], [2], [1], [2]]       
    end

    @testset "Quads" begin
        F = QuadFace{Int64}[[1,2,3,4],[3,4,5,6]]        
        con_E2F = con_edge_face(F)
        @test con_E2F == [[1], [1, 2], [1], [2], [2], [1], [2]]       
    end
end


@testset "con_face_face" verbose = true begin    
    @testset "Single triangle" begin
        F = TriangleFace{Int64}[[1,2,3]]
        E = meshedges(F;unique_only=false)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
        con_E2F = con_edge_face(F)
        con_F2E = con_face_edge(F)

        
        con_F2F = con_face_face(F)
        @test con_F2F == [[]]

        con_F2F_2 = con_face_face(F,E_uni,indReverse,con_E2F,con_F2E)
        @test con_F2F == con_F2F_2
    end

    @testset "Single quad" begin
        F = QuadFace{Int64}[[1,2,3,4]]
        E = meshedges(F;unique_only=false)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
        con_E2F = con_edge_face(F)
        con_F2E = con_face_edge(F)

        con_F2F = con_face_face(F)
        @test con_F2F == [[]]

        con_F2F_2 = con_face_face(F,E_uni,indReverse,con_E2F,con_F2E)
        @test con_F2F == con_F2F_2
    end

    @testset "Triangles" begin
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]        
        con_F2F = con_face_face(F)
        @test con_F2F == [[2], [1]]       
    end

    @testset "Quads" begin
        F = QuadFace{Int64}[[1,2,3,4],[3,4,5,6]]        
        con_F2F = con_face_face(F)
        @test con_F2F == [[2], [1]]       
    end
end


@testset "con_face_face_v" verbose = true begin    
    @testset "Single triangle" begin
        F = TriangleFace{Int64}[[1,2,3]]
        V = [Point3{Float64}(rand(3)) for _ in 1:3]        
        con_V2F = con_vertex_face(F,V) 

        con_F2F = con_face_face_v(F)
        @test con_F2F == [[]]

        con_F2F_2 = con_face_face_v(F,V,con_V2F)
        @test con_F2F == con_F2F_2
    end

    @testset "Single quad" begin
        F = QuadFace{Int64}[[1,2,3,4]]
        V = [Point3{Float64}(rand(4)) for _ in 1:4]        
        con_V2F = con_vertex_face(F,V)

        con_F2F = con_face_face_v(F)
        @test con_F2F == [[]]

        con_F2F_2 = con_face_face_v(F,V,con_V2F)
        @test con_F2F == con_F2F_2
    end

    @testset "Triangles" begin
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]        
        con_F2F = con_face_face_v(F)
        @test con_F2F == [[2], [1]]       
    end

    @testset "Quads" begin
        F = QuadFace{Int64}[[1,2,3,4],[3,4,5,6]]        
        con_F2F = con_face_face_v(F)
        @test con_F2F == [[2], [1]]       
    end
end


@testset "con_vertex_simplex" verbose = true begin    
    @testset "Single triangle" begin
        F = TriangleFace{Int64}[[1,2,3]]
        V = [Point3{Float64}(rand(3)) for _ in 1:3]        

        con_V2F = con_vertex_simplex(F)
        @test con_V2F == fill([1],length(F[1]))

        con_V2F_2 = con_vertex_simplex(F,V)
        @test con_V2F == con_V2F_2
    end

    @testset "Single quad" begin
        F = QuadFace{Int64}[[1,2,3,4]]
        V = [Point3{Float64}(rand(4)) for _ in 1:4]        

        con_V2F = con_vertex_simplex(F)
        @test con_V2F == fill([1],length(F[1]))

        con_V2F_2 = con_vertex_simplex(F,V)
        @test con_V2F == con_V2F_2
    end

    @testset "Triangles" begin
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]        
        con_V2F = con_vertex_simplex(F)
        @test con_V2F == [[1], [1, 2], [1, 2], [2]]      
    end

    @testset "Quads" begin
        F = QuadFace{Int64}[[1,2,3,4],[3,4,5,6]]        
        con_V2F = con_vertex_simplex(F)
        @test con_V2F == [[1], [1], [1, 2], [1, 2], [2], [2]]  
    end
end


@testset "con_vertex_face" verbose = true begin    
    @testset "Single triangle" begin
        F = TriangleFace{Int64}[[1,2,3]]
        V = [Point3{Float64}(rand(3)) for _ in 1:3]        

        con_V2F = con_vertex_face(F)
        @test con_V2F == fill([1],length(F[1]))

        con_V2F_2 = con_vertex_face(F,V)
        @test con_V2F == con_V2F_2
    end

    @testset "Single quad" begin
        F = QuadFace{Int64}[[1,2,3,4]]
        V = [Point3{Float64}(rand(4)) for _ in 1:4]        

        con_V2F = con_vertex_face(F)
        @test con_V2F == fill([1],length(F[1]))

        con_V2F_2 = con_vertex_face(F,V)
        @test con_V2F == con_V2F_2
    end

    @testset "Triangles" begin
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]        
        con_V2F = con_vertex_face(F)
        @test con_V2F == [[1], [1, 2], [1, 2], [2]]      
    end

    @testset "Quads" begin
        F = QuadFace{Int64}[[1,2,3,4],[3,4,5,6]]        
        con_V2F = con_vertex_face(F)
        @test con_V2F == [[1], [1], [1, 2], [1, 2], [2], [2]]  
    end
end


@testset "con_vertex_edge" verbose = true begin    
    @testset "Single triangle" begin
        F = TriangleFace{Int64}[[1,2,3]]
        E = meshedges(F)
        V = [Point3{Float64}(rand(3)) for _ in 1:3]        

        con_V2E = con_vertex_edge(E)
        @test con_V2E == [[1, 3], [1, 2], [2, 3]]

        con_V2E_2 = con_vertex_edge(E,V)
        @test con_V2E == con_V2E_2
    end

    @testset "Single quad" begin
        F = QuadFace{Int64}[[1,2,3,4]]
        E = meshedges(F)
        V = [Point3{Float64}(rand(4)) for _ in 1:4]        

        con_V2E = con_vertex_edge(E)
        @test con_V2E == [[1, 4], [1, 2], [2, 3], [3, 4]]

        con_V2E_2 = con_vertex_edge(E,V)
        @test con_V2E == con_V2E_2
    end

    @testset "Triangles" begin
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]     
        E = meshedges(F)  
        con_V2E = con_vertex_edge(E)
        @test con_V2E == [[1, 5], [1, 2, 3, 6], [2, 3, 4, 5], [4, 6]]    
    end

    @testset "Quads" begin
        F = QuadFace{Int64}[[1,2,3,4],[3,4,5,6]]        
        E = meshedges(F)
        con_V2E = con_vertex_edge(E)
        @test con_V2E == [[1, 7], [1, 3], [2, 3, 5, 8], [2, 4, 5, 7], [4, 6], [6, 8]]
    end
end


@testset "con_edge_edge" verbose = true begin    
    @testset "Single triangle" begin
        F = TriangleFace{Int64}[[1,2,3]]
        E = meshedges(F;unique_only=false)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
        con_V2E = con_vertex_edge(E_uni) 
                
        con_E2E = con_edge_edge(E_uni)
        @test con_E2E == [[3, 2], [1, 3], [2, 1]]

        con_E2E_2 = con_edge_edge(E_uni,con_V2E)
        @test con_E2E == con_E2E_2
    end

    @testset "Single quad" begin
        F = QuadFace{Int64}[[1,2,3,4]]
        E = meshedges(F;unique_only=false)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
        con_V2E = con_vertex_edge(E_uni) 

        con_E2E = con_edge_edge(E_uni)
        @test con_E2E == [[4, 2], [1, 3], [2, 4], [3, 1]]

        con_E2E_2 = con_edge_edge(E_uni,con_V2E)
        @test con_E2E == con_E2E_2
    end

    @testset "Triangles" begin
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]       
        E = meshedges(F;unique_only=false)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
        con_E2E = con_edge_edge(E_uni)
        @test con_E2E == [[4, 2, 5], [1, 5, 3, 4], [2, 4, 5], [2, 3, 1], [3, 1, 2]]      
    end

    @testset "Quads" begin
        F = QuadFace{Int64}[[1,2,3,4],[3,4,5,6]]      
        E = meshedges(F;unique_only=false)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)      
        con_E2E = con_edge_edge(E_uni)
        @test con_E2E == [[6, 3], [3, 7, 4, 6], [1, 2, 7], [2, 6, 5], [4, 7], [2, 4, 1], [5, 2, 3]]    
    end
end


@testset "con_vertex_vertex_f" verbose = true begin    
    @testset "Single triangle" begin
        F = TriangleFace{Int64}[[1,2,3]]        
        V = [Point3{Float64}(rand(3)) for _ in 1:3]        
        con_V2F = con_vertex_face(F)

        con_V2V = con_vertex_vertex_f(F)
        @test con_V2V == [[2, 3], [1, 3], [1, 2]]

        con_V2V_2 = con_vertex_vertex_f(F,V,con_V2F)
        @test con_V2V == con_V2V_2
    end

    @testset "Single quad" begin
        F = QuadFace{Int64}[[1,2,3,4]]        
        V = [Point3{Float64}(rand(4)) for _ in 1:4]        
        con_V2F = con_vertex_face(F)

        con_V2V = con_vertex_vertex_f(F)
        @test con_V2V == [[2, 3, 4], [1, 3, 4], [1, 2, 4], [1, 2, 3]]

        con_V2V_2 = con_vertex_vertex_f(F,V,con_V2F)
        @test con_V2V == con_V2V_2
    end

    @testset "Triangles" begin
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]             
        con_V2V = con_vertex_vertex_f(F)
        @test con_V2V == [[2, 3], [1, 3, 4], [1, 2, 4], [2, 3]]   
    end

    @testset "Quads" begin
        F = QuadFace{Int64}[[1,2,3,4],[3,4,5,6]]                
        con_V2V = con_vertex_vertex_f(F)
        @test con_V2V == [[2, 3, 4], [1, 3, 4], [1, 2, 4, 5, 6], [1, 2, 3, 5, 6], [3, 4, 6], [3, 4, 5]]
    end
end


@testset "con_vertex_vertex" verbose = true begin    
    @testset "Single triangle" begin
        F = TriangleFace{Int64}[[1,2,3]]        
        V = [Point3{Float64}(rand(3)) for _ in 1:3]        
        E = meshedges(F;unique_only=false)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
        con_V2E = con_vertex_edge(E_uni) 

        con_V2V = con_vertex_vertex(E)
        @test con_V2V == [[2, 3], [1, 3], [2, 1]]

        con_V2V_2 = con_vertex_vertex(E,V,con_V2E)
        @test con_V2V == con_V2V_2
    end

    @testset "Single quad" begin
        F = QuadFace{Int64}[[1,2,3,4]]        
        V = [Point3{Float64}(rand(4)) for _ in 1:4]        
        E = meshedges(F;unique_only=false)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)    
        con_V2E = con_vertex_edge(E_uni) 

        con_V2V = con_vertex_vertex(E)
        @test con_V2V == [[2, 4], [1, 3], [2, 4], [3, 1]]

        con_V2V_2 = con_vertex_vertex(E,V,con_V2E)
        @test con_V2V == con_V2V_2
    end

    @testset "Triangles" begin
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]      
        E = meshedges(F;unique_only=false)       
        con_V2V = con_vertex_vertex(E)
        @test con_V2V == [[2, 3], [1, 3, 4], [2, 4, 1], [3, 2]]
    end

    @testset "Quads" begin
        F = QuadFace{Int64}[[1,2,3,4],[3,4,5,6]]                
        E = meshedges(F;unique_only=false)
        con_V2V = con_vertex_vertex(E)
        @test con_V2V == [[2, 4], [1, 3], [4, 2, 6], [3, 5, 1], [4, 6], [5, 3]]
    end
end


@testset "meshconnectivity" verbose = true begin    
    @testset "Single triangle" begin
        F = TriangleFace{Int64}[[1,2,3]]        
        V = [Point3{Float64}(rand(3)) for _ in 1:3]        

        # EDGE-VERTEX connectivity
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)

        # FACE-EDGE connectivity
        con_F2E = con_face_edge(F,E_uni,indReverse)    

        # EDGE-FACE connectivity
        con_E2F = con_edge_face(F,E_uni)

        # FACE-FACE connectivity wrt edges
        con_F2F = con_face_face(F,E_uni,indReverse,con_E2F,con_F2E)

        # VERTEX-FACE connectivity
        con_V2F = con_vertex_face(F,V)

        # VERTEX-EDGE connectivity
        con_V2E = con_vertex_edge(E_uni,V)

        # EDGE-EDGE connectivity
        con_E2E = con_edge_edge(E_uni,con_V2E)

        # VERTEX-VERTEX connectivity wrt edges
        con_V2V = con_vertex_vertex(E_uni,V,con_V2E)

        # VERTEX-VERTEX connectivity wrt faces
        con_V2V_f = con_vertex_vertex_f(E_uni,V,con_V2E)

        # FACE-FACE connectivity wrt vertices
        con_F2F_v = con_face_face_v(F,con_V2F)
        
        C = meshconnectivity(F,V)

        @test typeof(C) == ConnectivitySet
        @test C.edge_vertex == E_uni
        @test C.edge_face == con_E2F
        @test C.edge_edge == con_E2E
        @test C.face_vertex == F
        @test C.face_edge == con_F2E
        @test C.face_face == con_F2F
        @test C.vertex_face == con_V2F
        @test C.vertex_vertex == con_V2V
        @test C.vertex_vertex_f == con_V2V_f
        @test C.face_face_v == con_F2F_v
    end

    @testset "Single quad" begin
        F = QuadFace{Int64}[[1,2,3,4]]        
        V = [Point3{Float64}(rand(4)) for _ in 1:4]   

        # EDGE-VERTEX connectivity
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)

        # FACE-EDGE connectivity
        con_F2E = con_face_edge(F,E_uni,indReverse)    

        # EDGE-FACE connectivity
        con_E2F = con_edge_face(F,E_uni)

        # FACE-FACE connectivity wrt edges
        con_F2F = con_face_face(F,E_uni,indReverse,con_E2F,con_F2E)

        # VERTEX-FACE connectivity
        con_V2F = con_vertex_face(F,V)

        # VERTEX-EDGE connectivity
        con_V2E = con_vertex_edge(E_uni,V)

        # EDGE-EDGE connectivity
        con_E2E = con_edge_edge(E_uni,con_V2E)

        # VERTEX-VERTEX connectivity wrt edges
        con_V2V = con_vertex_vertex(E_uni,V,con_V2E)

        # VERTEX-VERTEX connectivity wrt faces
        con_V2V_f = con_vertex_vertex_f(E_uni,V,con_V2E)

        # FACE-FACE connectivity wrt vertices
        con_F2F_v = con_face_face_v(F,con_V2F)

        C = meshconnectivity(F,V)

        @test typeof(C) == ConnectivitySet
        @test C.edge_vertex == E_uni
        @test C.edge_face == con_E2F
        @test C.edge_edge == con_E2E
        @test C.face_vertex == F
        @test C.face_edge == con_F2E
        @test C.face_face == con_F2F
        @test C.vertex_face == con_V2F
        @test C.vertex_vertex == con_V2V
        @test C.vertex_vertex_f == con_V2V_f
        @test C.face_face_v == con_F2F_v
    end

    @testset "Triangles" begin
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]      
        V = [Point3{Float64}(rand(4)) for _ in 1:4]   

        # EDGE-VERTEX connectivity
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)

        # FACE-EDGE connectivity
        con_F2E = con_face_edge(F,E_uni,indReverse)    

        # EDGE-FACE connectivity
        con_E2F = con_edge_face(F,E_uni)

        # FACE-FACE connectivity wrt edges
        con_F2F = con_face_face(F,E_uni,indReverse,con_E2F,con_F2E)

        # VERTEX-FACE connectivity
        con_V2F = con_vertex_face(F,V)

        # VERTEX-EDGE connectivity
        con_V2E = con_vertex_edge(E_uni,V)

        # EDGE-EDGE connectivity
        con_E2E = con_edge_edge(E_uni,con_V2E)

        # VERTEX-VERTEX connectivity wrt edges
        con_V2V = con_vertex_vertex(E_uni,V,con_V2E)

        # VERTEX-VERTEX connectivity wrt faces
        con_V2V_f = con_vertex_vertex_f(E_uni,V,con_V2E)

        # FACE-FACE connectivity wrt vertices
        con_F2F_v = con_face_face_v(F,con_V2F)

        C = meshconnectivity(F,V)

        @test typeof(C) == ConnectivitySet
        @test C.edge_vertex == E_uni
        @test C.edge_face == con_E2F
        @test C.edge_edge == con_E2E
        @test C.face_vertex == F
        @test C.face_edge == con_F2E
        @test C.face_face == con_F2F
        @test C.vertex_face == con_V2F
        @test C.vertex_vertex == con_V2V
        @test C.vertex_vertex_f == con_V2V_f
        @test C.face_face_v == con_F2F_v
    end

    @testset "Quads" begin
        F = QuadFace{Int64}[[1,2,3,4],[3,4,5,6]]                
        V = [Point3{Float64}(rand(4)) for _ in 1:6]   
        
        # EDGE-VERTEX connectivity
        E = meshedges(F)
        E_uni,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)

        # FACE-EDGE connectivity
        con_F2E = con_face_edge(F,E_uni,indReverse)    

        # EDGE-FACE connectivity
        con_E2F = con_edge_face(F,E_uni)

        # FACE-FACE connectivity wrt edges
        con_F2F = con_face_face(F,E_uni,indReverse,con_E2F,con_F2E)

        # VERTEX-FACE connectivity
        con_V2F = con_vertex_face(F,V)

        # VERTEX-EDGE connectivity
        con_V2E = con_vertex_edge(E_uni,V)

        # EDGE-EDGE connectivity
        con_E2E = con_edge_edge(E_uni,con_V2E)

        # VERTEX-VERTEX connectivity wrt edges
        con_V2V = con_vertex_vertex(E_uni,V,con_V2E)

        # VERTEX-VERTEX connectivity wrt faces
        con_V2V_f = con_vertex_vertex_f(E_uni,V,con_V2E)

        # FACE-FACE connectivity wrt vertices
        con_F2F_v = con_face_face_v(F,con_V2F)

        C = meshconnectivity(F,V)

        @test typeof(C) == ConnectivitySet
        @test C.edge_vertex == E_uni
        @test C.edge_face == con_E2F
        @test C.edge_edge == con_E2E
        @test C.face_vertex == F
        @test C.face_edge == con_F2E
        @test C.face_face == con_F2F
        @test C.vertex_face == con_V2F
        @test C.vertex_vertex == con_V2V
        @test C.vertex_vertex_f == con_V2V_f
        @test C.face_face_v == con_F2F_v
    end
end

@testset "mergevertices" begin
    eps_level = 1e-4
    r = sqrt(3)    
    M = cube(r)
    F = faces(M)
    V = coordinates(M)
    Fs,Vs = separate_vertices(F,V)
    Fm, Vm, indReverse = mergevertices(Fs, Vs; roundVertices = true)

    @test V isa Vector{Point3{Float64}}
    @test length(Vm) == length(V)
    @test isapprox(Vm[1], [-1.0, -1.0, -1.0], atol=eps_level)
    @test [indReverse[f] for f in Fm] == Fs # Reverse mapping 
    @test Fm isa Vector{QuadFace{Int64}} # Face type unaltered 
    @test length(Fm) == length(F) # Length is correct
    @test Fm[1] == [1, 2, 3, 4]

    Fm, Vm, indReverse = mergevertices(Fs, Vs; roundVertices = false)

    @test V isa Vector{Point3{Float64}}
    @test length(Vm) == length(V)
    @test isapprox(Vm[1], [-1.0, -1.0, -1.0], atol=eps_level)
    @test [indReverse[f] for f in Fm] == Fs # Reverse mapping 
    @test Fm isa Vector{QuadFace{Int64}} # Face type unaltered 
    @test length(Fm) == length(F) # Length is correct
    @test Fm[1] == [1, 2, 3, 4]
end


@testset "smoothmesh_laplacian" verbose = true begin

    eps_level = 1e-4
    M = tetrahedron(1.0)
    F = faces(M)
    V = coordinates(M)
    F,V = subtri(F,V,3)

    ind = round.(Int64,range(1,length(V),5))

    @testset "errors" begin
        λ = -0.5 # Laplacian smoothing parameter
        n = 10 # Number of iterations 
        @test_throws ArgumentError smoothmesh_laplacian(F,V,n,λ)

        λ = 0.5 # Laplacian smoothing parameter
        n = -3 # Number of iterations 
        @test_throws ArgumentError smoothmesh_laplacian(F,V,n,λ)
    end

    @testset "n=0" begin
        λ = 0.5 # Laplacian smoothing parameter
        n = 0 # Number of iterations 
        @test V == smoothmesh_laplacian(F,V,n,λ)
    end

    @testset "Unconstrained smoothing" begin
        λ = 0.5 # Laplacian smoothing parameter
        n = 10 # Number of iterations 
        Vs = smoothmesh_laplacian(F,V,n,λ)
    
        @test length(V) == length(Vs)
        @test isapprox(Vs[ind],Point3{Float64}[[-0.5267934833030736, -0.3041443593923703, -0.21506253898598338],
         [0.443754269985866, -0.2562016472303856, -0.07946007147312682], 
         [-0.08493319468029126, 0.11913157297777958, 0.4220011641990833], 
         [-6.987080667128806e-19, 0.1939578299475715, 0.3708898233176402], 
         [-0.022334842065796563, -0.012895027078995421, 0.6087149725933829]],
                                                atol=eps_level)
    end

    @testset "Constrained smoothing" begin
        λ = 0.5 # Laplacian smoothing parameter
        n = 10 # Number of iterations 
        Vs = smoothmesh_laplacian(F,V,n,λ; constrained_points = [5])
    
        @test Vs[5] == V[5]
        @test length(V) == length(Vs)
        @test isapprox(Vs[ind],Point3{Float64}[[-0.5267934833030736, -0.3045704088157244, -0.21536380142235773],
        [0.443754269985866, -0.2568254964220775, -0.07990119946700998], 
        [-0.08493319468029126, 0.11913135061947012, 0.4220010069680148], 
        [-6.987080667128806e-19, 0.19395781755288108, 0.3708898145532705], 
        [-0.022334842065796563, -0.012895869445908863, 0.6087143769500261]], atol=eps_level)
    end

end


@testset "smoothmesh_hc" verbose = true begin
    eps_level = 1e-4
    M = tetrahedron(1.0)
    F = faces(M)
    V = coordinates(M)
    F,V = subtri(F,V,3)

    ind = round.(Int64,range(1,length(V),5))

    @testset "errors" begin
        n=10
        α=-0.1
        β=0.5
        @test_throws ArgumentError smoothmesh_hc(F,V,n,α,β)

        n=10
        α=0.1
        β=-0.5
        @test_throws ArgumentError smoothmesh_hc(F,V,n,α,β)

        n=-3
        α=0.1
        β=0.5
        @test_throws ArgumentError smoothmesh_hc(F,V,n,α,β)
    end

    @testset "n=0" begin
        n=0
        α=0.1
        β=0.5
        @test V == smoothmesh_hc(F,V,n,α,β)
    end

    @testset "Unconstrained smoothing" begin
        n=10
        α=0.1
        β=0.5
        Vs = smoothmesh_hc(F,V,n,α,β)
    
        @test length(V) == length(Vs)
        @test isapprox(Vs[ind],Point3{Float64}[[-0.717172128187643, -0.4140595212644323, -0.2927842953009355], 
        [0.5758586318390522, -0.3324721361074463, -0.0620203220858017], 
        [-0.11011137811964947, 0.17615457164191245, 0.521908672591213], 
        [-8.683563963552147e-19, 0.2903204379104125, 0.4677777851758052], 
        [-0.04381485137375324, -0.025296516235139878, 0.8068035332217547]],atol=eps_level)
    end

    @testset "Constrained smoothing" begin
        n=10
        α=0.1
        β=0.5
        Vs = smoothmesh_hc(F,V,n,α,β; constrained_points = [5])
    
        @test Vs[5] == V[5]
    
        @test length(V) == length(Vs)
        @test isapprox(Vs[ind],Point3{Float64}[[-0.717172128187643, -0.4136411138967243, -0.29248843661393087], 
        [0.5758586318390522, -0.33191585522870704, -0.061626972104200775], 
        [-0.11011137811964947, 0.17615378734209927, 0.5219081180074967], 
        [-8.683563963552147e-19, 0.2903235073180716, 0.46777995557477525], 
        [-0.04381485137375324, -0.025308290081523535, 0.8067952078551363]],atol=eps_level)
    end

    @testset "Distance threshold termination" begin
        n=1000
        α=0.1
        β=0.5
        tolDist=1e-2
        Vs = smoothmesh_hc(F,V,n,α,β;  tolDist=tolDist)
    
        @test length(V) == length(Vs)
        @test isapprox(Vs[ind],Point3{Float64}[[-0.7077940615748582, -0.40864509198106286, -0.28895571563840994], 
        [0.5726730482498883, -0.3306329385647165, -0.06830402485575475], 
        [-0.10900177329770618, 0.17092653723203155, 0.5232535614837974], 
        [-2.4395652898844805e-18, 0.2818887131438471, 0.4682163057073023], 
        [-0.04159517654997694, -0.024014986378119194, 0.7989425080429475]],atol=eps_level)
    end
    
end


@testset "quadplate" begin
    eps_level = 1e-4

    plateDim = [5,5]
    plateElem = [4,3]
    F,V = quadplate(plateDim,plateElem)

    @test V isa Vector{Point3{Float64}}
    @test length(V) == prod(plateElem.+1)
    @test isapprox(V, Point3{Float64}[[-2.5, -2.5, 0.0], [-1.25, -2.5, 0.0], [0.0, -2.5, 0.0], [1.25, -2.5, 0.0], [2.5, -2.5, 0.0], [-2.5, -0.8333333333333334, 0.0], [-1.25, -0.8333333333333334, 0.0], [0.0, -0.8333333333333334, 0.0], [1.25, -0.8333333333333334, 0.0], [2.5, -0.8333333333333334, 0.0], [-2.5, 0.8333333333333334, 0.0], [-1.25, 0.8333333333333334, 0.0], [0.0, 0.8333333333333334, 0.0], [1.25, 0.8333333333333334, 0.0], [2.5, 0.8333333333333334, 0.0], [-2.5, 2.5, 0.0], [-1.25, 2.5, 0.0], [0.0, 2.5, 0.0], [1.25, 2.5, 0.0], [2.5, 2.5, 0.0]], atol=eps_level)
    @test F isa Vector{QuadFace{Int64}}
    @test length(F) == prod(plateElem)
    @test F[1] == [1, 2, 7, 6]
end


@testset "quadsphere" begin
    eps_level = 1e-4

    r = 1.0
    F, V = quadsphere(3, r)
    
    @test V isa Vector{Point3{Float64}}
    @test length(V) == 386
    @test isapprox(V[1], [-0.5773502691896258, -0.5773502691896258, -0.5773502691896258], atol=eps_level)
    @test F isa Vector{QuadFace{Int64}}
    @test length(F) == 384
    @test F[1] == [1, 99, 291, 187]
end

@testset "simplex2vertexdata" verbose = true begin
    eps_level = 1e-4

    # Single face/element
    F1 = [[1,2,3,4,5,6]]
    V1 = [GeometryBasics.Point3(rand(3)) for _=1:length(F1[1])]

    # A quad mesh featuring a variation in terms of face areas and vertex connectivity 
    Mq = cube(1.0)
    Fq = faces(Mq)
    Vq = coordinates(Mq)
    Fq,Vq = subquad(Fq,Vq,1; method=:Catmull_Clark)

    # A triangle mesh featuring a variation in terms of face areas and vertex connectivity 
    Mt = tetrahedron(1.0)
    Ft = faces(Mt)
    Vt = coordinates(Mt)
    Ft,Vt = subtri(Ft,Vt,1; method=:Loop)

    # A hexahedral mesh 
    Eh,Vh,_,_,_ = hexbox([1.0,1.0,1.0],[2,2,2])

    @testset "Errors" begin
        DF = [1.0] # Face data (here face number)
        @test_throws ArgumentError simplex2vertexdata(F1,DF,nothing; weighting=:area)
    end

    @testset "Vector Float64 data, weighting=:none" begin
        # Single element
        DF = [1.0] # Face data (here face number)
        DV = simplex2vertexdata(F1,DF,V1; weighting=:none)
        @test isapprox(DV,ones(Float64,length(F1[1])), atol=eps_level)

        # Quads
        DF = collect(Float64,1:length(Fq)) # Face data (here face numbers)
        DV = simplex2vertexdata(Fq,DF,Vq; weighting=:none)
        @test isapprox(DV,[12.0, 9.666666666666666, 12.666666666666666, 15.666666666666666, 13.0, 10.0, 12.333333333333334, 14.666666666666666, 6.5, 11.5, 8.5, 10.0, 14.0, 9.0, 12.5, 16.5, 20.5, 16.5, 11.5, 13.0, 2.5, 6.5, 10.5, 14.5, 18.5, 22.5], atol=eps_level)

        # Triangles
        DF = collect(Float64,1:length(Ft)) # Face data (here face numbers)
        DV = simplex2vertexdata(Ft,DF,Vt; weighting=:none)
        @test isapprox(DV,[10.333333333333334, 11.333333333333334, 13.333333333333334, 7.0, 6.833333333333333, 7.166666666666667, 8.166666666666666, 7.666666666666667, 8.666666666666666, 8.5], atol=eps_level)

        # Hexahedral elements
        DF = collect(Float64,1:length(Eh)) # Face data (here face numbers)
        DV = simplex2vertexdata(Eh,DF,Vh; weighting=:none)
        @test isapprox(DV,[1.0, 1.5, 2.0, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0, 3.0, 3.5, 4.0, 4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.0, 5.5, 6.0, 6.0, 6.5, 7.0, 7.0, 7.5, 8.0], atol=eps_level)
    end
    
    @testset "Vector Float64 data, weighting=:area" begin
        # Single element
        DF = [1.0] # Face data (here face number)
        DV = simplex2vertexdata(F1,DF,V1; weighting=:area)
        @test isapprox(DV,ones(Float64,length(F1[1])), atol=eps_level)

        # Quads
        DF = collect(Float64,1:length(Fq)) # Face data (here face numbers)
        DV = simplex2vertexdata(Fq,DF,Vq; weighting=:area)
        @test isapprox(DV,[12.0, 9.666666666666668, 12.666666666666666, 15.666666666666664, 13.0, 10.0, 12.333333333333332, 14.666666666666666, 6.500000000000001, 11.5, 8.5, 10.0, 14.0, 9.0, 12.5, 16.5, 20.5, 16.5, 11.5, 13.000000000000002, 2.5, 6.500000000000001, 10.5, 14.5, 18.5, 22.5], atol=eps_level)

        # Triangles
        DF = collect(Float64,1:length(Ft)) # Face data (here face numbers)
        DV = simplex2vertexdata(Ft,DF,Vt; weighting=:area)
        @test isapprox(DV,[10.333333333333334, 11.333333333333334, 13.333333333333336, 6.999999999999999, 5.095917942265425, 5.646428199482246, 6.646428199482247, 6.146428199482247, 6.49489742783178, 6.545407685048602], atol=eps_level)

    end
    
    @testset "Vector of Vector Float64 data, weighting=:none" begin

        # Single element
        DF = [[1.0,2.0,π]] # Face data (here face number)
        DV = simplex2vertexdata(F1,DF,V1; weighting=:none)
        @test isapprox(DV,repeat(DF,length(V1)), atol=eps_level)

        # Quads
        DF = [i.*[1.0,2.0,3.0] for i in eachindex(Fq)] # Vector data for each face
        DV = simplex2vertexdata(Fq,DF,Vq; weighting=:none)
        @test isapprox(DV,[[12.0, 24.0, 36.0], [9.666666666666666, 19.333333333333332, 29.0], [12.666666666666666, 25.333333333333332, 38.0], [15.666666666666666, 31.333333333333332, 47.0], [13.0, 26.0, 39.0], [10.0, 20.0, 30.0], [12.333333333333334, 24.666666666666668, 37.0], [14.666666666666666, 29.333333333333332, 44.0], [6.5, 13.0, 19.5], [11.5, 23.0, 34.5], [8.5, 17.0, 25.5], [10.0, 20.0, 30.0], [14.0, 28.0, 42.0], [9.0, 18.0, 27.0], [12.5, 25.0, 37.5], [16.5, 33.0, 49.5], [20.5, 41.0, 61.5], [16.5, 33.0, 49.5], [11.5, 23.0, 34.5], [13.0, 26.0, 39.0], [2.5, 5.0, 7.5], [6.5, 13.0, 19.5], [10.5, 21.0, 31.5], [14.5, 29.0, 43.5], [18.5, 37.0, 55.5], [22.5, 45.0, 67.5]], atol=eps_level)

        # Triangles
        DF = [i.*[1.0,2.0,3.0] for i in eachindex(Ft)] # Vector data for each face
        DV = simplex2vertexdata(Ft,DF,Vt; weighting=:none)
        @test isapprox(DV,[[10.333333333333334, 20.666666666666668, 31.0], [11.333333333333334, 22.666666666666668, 34.0], [13.333333333333334, 26.666666666666668, 40.0], [7.0, 14.0, 21.0], [6.833333333333333, 13.666666666666666, 20.5], [7.166666666666667, 14.333333333333334, 21.5], [8.166666666666666, 16.333333333333332, 24.5], [7.666666666666667, 15.333333333333334, 23.0], [8.666666666666666, 17.333333333333332, 26.0], [8.5, 17.0, 25.5]], atol=eps_level)

        # Hexahedral elements
        DF = [i.*[1.0,2.0,3.0] for i in eachindex(Eh)] # Vector data for each face
        DV = simplex2vertexdata(Eh,DF,Vh; weighting=:none)
        @test isapprox(DV,[[1.0, 2.0, 3.0], [1.5, 3.0, 4.5], [2.0, 4.0, 6.0], [2.0, 4.0, 6.0], [2.5, 5.0, 7.5], [3.0, 6.0, 9.0], [3.0, 6.0, 9.0], [3.5, 7.0, 10.5], [4.0, 8.0, 12.0], [3.0, 6.0, 9.0], [3.5, 7.0, 10.5], [4.0, 8.0, 12.0], [4.0, 8.0, 12.0], [4.5, 9.0, 13.5], [5.0, 10.0, 15.0], [5.0, 10.0, 15.0], [5.5, 11.0, 16.5], [6.0, 12.0, 18.0], [5.0, 10.0, 15.0], [5.5, 11.0, 16.5], [6.0, 12.0, 18.0], [6.0, 12.0, 18.0], [6.5, 13.0, 19.5], [7.0, 14.0, 21.0], [7.0, 14.0, 21.0], [7.5, 15.0, 22.5], [8.0, 16.0, 24.0]], atol=eps_level)
    end

    @testset "Vector of Matrix Float64 data, weighting=:none" begin

        # Single element
        DF = [[1.0 2.0 3.0; 4.0 5.0 6.0]] # Face data (here face number)
        DV = simplex2vertexdata(F1,DF,V1; weighting=:none)
        @test isapprox(DV,repeat(DF,length(V1)), atol=eps_level)

        # Quads
        DF = [i.*[1.0 2.0 3.0; 4.0 5.0 6.0] for i in eachindex(Fq)] # Matrix data for each face
        DV = simplex2vertexdata(Fq,DF,Vq; weighting=:area)
        @test isapprox(DV,[[12.0 24.0 36.0; 48.0 60.0 72.0], [9.666666666666668 19.333333333333336 28.999999999999996; 38.66666666666667 48.33333333333333 57.99999999999999], [12.666666666666666 25.333333333333332 38.0; 50.666666666666664 63.33333333333333 76.0], [15.666666666666664 31.33333333333333 47.0; 62.66666666666666 78.33333333333333 94.0], [13.0 26.0 39.0; 52.0 65.0 78.0], [10.0 20.0 30.0; 40.0 50.0 60.0], [12.333333333333332 24.666666666666664 37.0; 49.33333333333333 61.666666666666664 74.0], [14.666666666666666 29.333333333333332 44.0; 58.666666666666664 73.33333333333334 88.0], [6.500000000000001 13.000000000000002 19.5; 26.000000000000004 32.5 39.0], [11.5 23.0 34.5; 46.0 57.5 69.0], [8.5 17.0 25.500000000000004; 34.0 42.5 51.00000000000001], [10.0 20.0 30.0; 40.0 50.0 60.0], [14.0 28.0 42.0; 56.0 70.0 84.0], [9.0 18.0 27.0; 36.0 44.99999999999999 54.0], [12.5 25.0 37.50000000000001; 50.0 62.5 75.00000000000001], [16.5 33.0 49.5; 66.0 82.5 99.0], [20.5 41.0 61.5; 82.0 102.50000000000001 123.0], [16.5 33.0 49.5; 66.0 82.5 99.0], [11.5 23.0 34.5; 46.0 57.5 69.0], [13.000000000000002 26.000000000000004 39.0; 52.00000000000001 65.0 78.0], [2.5 5.0 7.5; 10.0 12.5 15.0], [6.500000000000001 13.000000000000002 19.5; 26.000000000000004 32.5 39.0], [10.5 21.0 31.5; 42.0 52.49999999999999 63.0], [14.5 29.0 43.5; 58.0 72.49999999999999 87.0], [18.5 37.0 55.5; 74.0 92.5 111.0], [22.5 45.0 67.5; 90.0 112.50000000000001 135.0]], atol=eps_level)

        # Triangles
        DF = [i.*[1.0 2.0 3.0; 4.0 5.0 6.0] for i in eachindex(Ft)] # Matrix data for each face
        DV = simplex2vertexdata(Ft,DF,Vt; weighting=:none)
        @test isapprox(DV,[[10.333333333333334 20.666666666666668 31.0; 41.333333333333336 51.666666666666664 62.0], [11.333333333333334 22.666666666666668 34.0; 45.333333333333336 56.666666666666664 68.0], [13.333333333333334 26.666666666666668 40.0; 53.333333333333336 66.66666666666667 80.0], [7.0 14.0 21.0; 28.0 35.0 42.0], [6.833333333333333 13.666666666666666 20.5; 27.333333333333332 34.166666666666664 41.0], [7.166666666666667 14.333333333333334 21.5; 28.666666666666668 35.833333333333336 43.0], [8.166666666666666 16.333333333333332 24.5; 32.666666666666664 40.833333333333336 49.0], [7.666666666666667 15.333333333333334 23.0; 30.666666666666668 38.333333333333336 46.0], [8.666666666666666 17.333333333333332 26.0; 34.666666666666664 43.333333333333336 52.0], [8.5 17.0 25.5; 34.0 42.5 51.0]], atol=eps_level)

        # Hexahedral elements
        DF = [i.*[1.0 2.0 3.0; 4.0 5.0 6.0] for i in eachindex(Eh)] # Matrix data for each face
        DV = simplex2vertexdata(Eh,DF,Vh; weighting=:none)
        @test isapprox(DV,[[1.0 2.0 3.0; 4.0 5.0 6.0], [1.5 3.0 4.5; 6.0 7.5 9.0], [2.0 4.0 6.0; 8.0 10.0 12.0], [2.0 4.0 6.0; 8.0 10.0 12.0], [2.5 5.0 7.5; 10.0 12.5 15.0], [3.0 6.0 9.0; 12.0 15.0 18.0], [3.0 6.0 9.0; 12.0 15.0 18.0], [3.5 7.0 10.5; 14.0 17.5 21.0], [4.0 8.0 12.0; 16.0 20.0 24.0], [3.0 6.0 9.0; 12.0 15.0 18.0], [3.5 7.0 10.5; 14.0 17.5 21.0], [4.0 8.0 12.0; 16.0 20.0 24.0], [4.0 8.0 12.0; 16.0 20.0 24.0], [4.5 9.0 13.5; 18.0 22.5 27.0], [5.0 10.0 15.0; 20.0 25.0 30.0], [5.0 10.0 15.0; 20.0 25.0 30.0], [5.5 11.0 16.5; 22.0 27.5 33.0], [6.0 12.0 18.0; 24.0 30.0 36.0], [5.0 10.0 15.0; 20.0 25.0 30.0], [5.5 11.0 16.5; 22.0 27.5 33.0], [6.0 12.0 18.0; 24.0 30.0 36.0], [6.0 12.0 18.0; 24.0 30.0 36.0], [6.5 13.0 19.5; 26.0 32.5 39.0], [7.0 14.0 21.0; 28.0 35.0 42.0], [7.0 14.0 21.0; 28.0 35.0 42.0], [7.5 15.0 22.5; 30.0 37.5 45.0], [8.0 16.0 24.0; 32.0 40.0 48.0]], atol=eps_level)
    end    
    
end


@testset "vertex2simplexdata" verbose = true begin
    eps_level = 1e-4

    # Single face/element
    F1 = [[1,2,3,4,5,6]]
    V1 = [GeometryBasics.Point3(rand(3)) for _=1:length(F1[1])]

    # A quad mesh featuring a variation in terms of face areas and vertex connectivity 
    Mq = cube(1.0)
    Fq = faces(Mq)
    Vq = coordinates(Mq)
    Fq,Vq = subquad(Fq,Vq,1; method=:Catmull_Clark)

    # A triangle mesh featuring a variation in terms of face areas and vertex connectivity 
    Mt = tetrahedron(1.0)
    Ft = faces(Mt)
    Vt = coordinates(Mt)
    Ft,Vt = subtri(Ft,Vt,1; method=:Loop)

    # A hexahedral mesh 
    Eh,Vh,_,_,_ = hexbox([1.0,1.0,1.0],[2,2,2])

    @testset "Vector Float64 data" begin
        # Single element
        DV = collect(Float64,1:length(V1)) # Vertex data (here node number)
        DF = vertex2simplexdata(F1,DV)
        @test isapprox(DF,[mean(DV)], atol=eps_level)

        # Quads
        DV = collect(Float64,1:length(Vq)) # Vertex data (here node number)
        DF = vertex2simplexdata(Fq,DV)
        @test isapprox(DF,[12.75, 11.5, 14.25, 16.0, 13.25, 12.75, 12.75, 12.75, 14.25, 13.75, 12.25, 12.75, 14.25, 14.75, 14.25, 13.75, 14.5, 15.0, 16.25, 15.75, 16.0, 15.5, 16.25, 16.75], atol=eps_level)

        # Triangles
        DV = collect(Float64,1:length(Vt)) # Vertex data (here node number)
        DF = vertex2simplexdata(Ft,DV)
        @test isapprox(DF,[8.0, 6.333333333333333, 7.333333333333333, 8.333333333333334, 5.333333333333333, 6.0, 5.666666666666667, 6.333333333333333, 5.333333333333333, 4.333333333333333, 6.333333333333333, 6.333333333333333, 7.333333333333333, 4.666666666666667, 5.666666666666667, 6.666666666666667], atol=eps_level)

        # Hexahedral elements
        DV = collect(Float64,1:length(Vh)) # Vertex data (here node number)
        DF = vertex2simplexdata(Eh,DV)
        @test isapprox(DF,[7.5, 8.5, 10.5, 11.5, 16.5, 17.5, 19.5, 20.5], atol=eps_level)
    end
    
    @testset "Vector of Vector Float64 data" begin

        # Single element
        DV = [i*[1.0,2.0,π] for i in eachindex(V1)] # Vector data for each vertex
        DF = vertex2simplexdata(F1,DV)
        @test isapprox(DF,[mean(DV)], atol=eps_level)

        # Quads
        DV = [i.*[1.0,2.0,π] for i in eachindex(Vq)] # Vector data for each vertex
        DF = vertex2simplexdata(Fq,DV)
        @testset "simplexcenter" begin
            eps_level = 1e-4
            F, V = geosphere(2, 1.0)
            VC = simplexcenter(F, V)
        
            @test VC isa typeof(V)
            @test length(VC) == length(F)
            @test isapprox(VC[1:30:end], Point3{Float64}[[-0.3504874080794224, 0.0, -0.9175879469807824], [0.9252211181650858, 0.17425248968910703, -0.2898716471939399], [-0.17425248968910703, -0.2898716471939399, -0.9252211181650858], [0.870241439674047, 0.4295322227262335, -0.15715894749713352], [-0.0876218520198556, 0.14524212567637496, -0.9709982913596904], [0.9709982913596904, 0.0876218520198556, 0.14524212567637496], [-0.5759258522984322, -0.1565463433383235, -0.7844827600958122], [0.7844827600958122, 0.5759258522984322, -0.1565463433383235], [-0.4295322227262335, -0.15715894749713352, -0.870241439674047], [0.8917525488507145, 0.08663063766925146, -0.41096206852816675], [-0.08663063766925146, -0.41096206852816675, -0.8917525488507145]], atol=eps_level)
        end
        @test isapprox(DF,[[12.75, 25.5, 40.05530633326987], [11.5, 23.0, 36.12831551628262], [14.25, 28.5, 44.76769531365455], [16.0, 32.0, 50.26548245743669], [13.25, 26.5, 41.62610266006476], [12.75, 25.5, 40.05530633326986], [12.75, 25.5, 40.055306333269854], [12.75, 25.5, 40.05530633326986], [14.25, 28.5, 44.76769531365455], [13.75, 27.5, 43.19689898685965], [12.25, 24.5, 38.48451000647496], [12.75, 25.5, 40.05530633326986], [14.25, 28.5, 44.767695313654556], [14.75, 29.5, 46.33849164044945], [14.25, 28.5, 44.767695313654556], [13.75, 27.5, 43.19689898685965], [14.5, 29.0, 45.553093477052], [15.0, 30.0, 47.1238898038469], [16.25, 32.5, 51.050880620834135], [15.75, 31.5, 49.48008429403924], [16.0, 32.0, 50.26548245743669], [15.5, 31.0, 48.69468613064179], [16.25, 32.5, 51.05088062083414], [16.75, 33.5, 52.62167694762904]], atol=eps_level)

        # Triangles
        DV = [i.*[1.0,2.0,π] for i in eachindex(Vt)] # Vector data for each vertex
        DF = vertex2simplexdata(Ft,DV)
        @test isapprox(DF,[[8.0, 16.0, 25.132741228718345], [6.333333333333333, 12.666666666666666, 19.896753472735355], [7.333333333333333, 14.666666666666666, 23.03834612632515], [8.333333333333334, 16.666666666666668, 26.179938779914945], [5.333333333333333, 10.666666666666666, 16.755160819145562], [6.0, 12.0, 18.84955592153876], [5.666666666666667, 11.333333333333334, 17.802358370342162], [6.333333333333333, 12.666666666666666, 19.896753472735355], [5.333333333333333, 10.666666666666666, 16.755160819145562], [4.333333333333333, 8.666666666666666, 13.613568165555769], [6.333333333333333, 12.666666666666666, 19.896753472735355], [6.333333333333333, 12.666666666666666, 19.896753472735355], [7.333333333333333, 14.666666666666666, 23.03834612632515], [4.666666666666667, 9.333333333333334, 14.66076571675237], [5.666666666666667, 11.333333333333334, 17.80235837034216], [6.666666666666667, 13.333333333333334, 20.943951023931955]], atol=eps_level)

        # Hexahedral elements
        DV = [i.*[1.0,2.0,π] for i in eachindex(Vh)] # Vector data for each face
        DF = vertex2simplexdata(Ft,DV)
        @test isapprox(DF,[[8.0, 16.0, 25.132741228718345], [6.333333333333333, 12.666666666666666, 19.896753472735355], [7.333333333333333, 14.666666666666666, 23.03834612632515], [8.333333333333334, 16.666666666666668, 26.179938779914945], [5.333333333333333, 10.666666666666666, 16.755160819145562], [6.0, 12.0, 18.84955592153876], [5.666666666666667, 11.333333333333334, 17.802358370342162], [6.333333333333333, 12.666666666666666, 19.896753472735355], [5.333333333333333, 10.666666666666666, 16.755160819145562], [4.333333333333333, 8.666666666666666, 13.613568165555769], [6.333333333333333, 12.666666666666666, 19.896753472735355], [6.333333333333333, 12.666666666666666, 19.896753472735355], [7.333333333333333, 14.666666666666666, 23.03834612632515], [4.666666666666667, 9.333333333333334, 14.66076571675237], [5.666666666666667, 11.333333333333334, 17.80235837034216], [6.666666666666667, 13.333333333333334, 20.943951023931955]], atol=eps_level)
    end

    @testset "Vector of Matrix Float64 data" begin

        # Single element
        DV = [i*[1.0 2.0 3.0; 4.0 5.0 6.0] for i in eachindex(V1)] # Matrix data for each vertex
        DF = vertex2simplexdata(F1,DV)
        @test isapprox(DF,[mean(DV)], atol=eps_level)

        # Quads
        DV = [i.*[1.0 2.0 3.0; 4.0 5.0 6.0] for i in eachindex(Vq)] # Matrix data for each vertex
        DF = vertex2simplexdata(Fq,DV)
        @test isapprox(DF,[[12.75 25.5 38.25; 51.0 63.75 76.5], [11.5 23.0 34.5; 46.0 57.5 69.0], [14.25 28.5 42.75; 57.0 71.25 85.5], [16.0 32.0 48.0; 64.0 80.0 96.0], [13.25 26.5 39.75; 53.0 66.25 79.5], [12.75 25.5 38.25; 51.0 63.75 76.5], [12.75 25.5 38.25; 51.0 63.75 76.5], [12.75 25.5 38.25; 51.0 63.75 76.5], [14.25 28.5 42.75; 57.0 71.25 85.5], [13.75 27.5 41.25; 55.0 68.75 82.5], [12.25 24.5 36.75; 49.0 61.25 73.5], [12.75 25.5 38.25; 51.0 63.75 76.5], [14.25 28.5 42.75; 57.0 71.25 85.5], [14.75 29.5 44.25; 59.0 73.75 88.5], [14.25 28.5 42.75; 57.0 71.25 85.5], [13.75 27.5 41.25; 55.0 68.75 82.5], [14.5 29.0 43.5; 58.0 72.5 87.0], [15.0 30.0 45.0; 60.0 75.0 90.0], [16.25 32.5 48.75; 65.0 81.25 97.5], [15.75 31.5 47.25; 63.0 78.75 94.5], [16.0 32.0 48.0; 64.0 80.0 96.0], [15.5 31.0 46.5; 62.0 77.5 93.0], [16.25 32.5 48.75; 65.0 81.25 97.5], [16.75 33.5 50.25; 67.0 83.75 100.5]], atol=eps_level)

        # Triangles
        DV = [i.*[1.0 2.0 3.0; 4.0 5.0 6.0] for i in eachindex(Vt)] # Matrix data for each vertex
        DF = vertex2simplexdata(Ft,DV)
        @test isapprox(DF,[[8.0 16.0 24.0; 32.0 40.0 48.0], [6.333333333333333 12.666666666666666 19.0; 25.333333333333332 31.666666666666668 38.0], [7.333333333333333 14.666666666666666 22.0; 29.333333333333332 36.666666666666664 44.0], [8.333333333333334 16.666666666666668 25.0; 33.333333333333336 41.666666666666664 50.0], [5.333333333333333 10.666666666666666 16.0; 21.333333333333332 26.666666666666668 32.0], [6.0 12.0 18.0; 24.0 30.0 36.0], [5.666666666666667 11.333333333333334 17.0; 22.666666666666668 28.333333333333332 34.0], [6.333333333333333 12.666666666666666 19.0; 25.333333333333332 31.666666666666668 38.0], [5.333333333333333 10.666666666666666 16.0; 21.333333333333332 26.666666666666668 32.0], [4.333333333333333 8.666666666666666 13.0; 17.333333333333332 21.666666666666668 26.0], [6.333333333333333 12.666666666666666 19.0; 25.333333333333332 31.666666666666668 38.0], [6.333333333333333 12.666666666666666 19.0; 25.333333333333332 31.666666666666668 38.0], [7.333333333333333 14.666666666666666 22.0; 29.333333333333332 36.666666666666664 44.0], [4.666666666666667 9.333333333333334 14.0; 18.666666666666668 23.333333333333332 28.0], [5.666666666666667 11.333333333333334 17.0; 22.666666666666668 28.333333333333332 34.0], [6.666666666666667 13.333333333333334 20.0; 26.666666666666668 33.333333333333336 40.0]], atol=eps_level)

        # Hexahedral elements
        DV = [i.*[1.0 2.0 3.0; 4.0 5.0 6.0] for i in eachindex(Vh)] # Matrix data for each face
        DF = vertex2simplexdata(Eh,DV)
        @test isapprox(DF,[[7.5 15.0 22.5; 30.0 37.5 45.0], [8.5 17.0 25.5; 34.0 42.5 51.0], [10.5 21.0 31.5; 42.0 52.5 63.0], [11.5 23.0 34.5; 46.0 57.5 69.0], [16.5 33.0 49.5; 66.0 82.5 99.0], [17.5 35.0 52.5; 70.0 87.5 105.0], [19.5 39.0 58.5; 78.0 97.5 117.0], [20.5 41.0 61.5; 82.0 102.5 123.0]], atol=eps_level)
    end        
end


@testset "simplexcenter" verbose = true  begin
    eps_level = 1e-4

    @testset "Triangles" begin
        F, V = geosphere(2, 1.0)
        VC = simplexcenter(F, V)
        ind = round.(Int64,range(1,length(VC),5))

        @test VC isa typeof(V)
        @test length(VC) == length(F)
        @test isapprox(VC[ind], Point3{Float64}[[-0.3504874080794224, 0.0, -0.9175879469807824], [-0.4295322227262335, 0.15715894749713352, -0.870241439674047], [-0.7866254783422401, 0.3950724390612581, -0.4435636037400384], [-0.6300791350039167, 0.29832147806366355, -0.6968609080759566], [-0.805121911181463, 0.14017131621592582, -0.5511333847440926]], atol=eps_level)
    end

    @testset "Quadrilaterals" begin
        F, V = quadsphere(2, 1.0)
        VC = simplexcenter(F, V)
        ind = round.(Int64,range(1,length(VC),5))
        @test VC isa typeof(V)
        @test length(VC) == length(F)
        @test isapprox(VC[ind], Point3{Float64}[[-0.4802860138667546, -0.4802860138667546, -0.6949720954766154], [-0.4802860138667546, 0.4802860138667546, 0.6949720954766154], [-0.7899092719339054, -0.16747661189958585, -0.5326696383253359], [0.7899092719339054, -0.16747661189958585, 0.5326696383253359], [0.16747661189958585, -0.7899092719339054, -0.5326696383253359]], atol=eps_level)
    end
end


@testset "normalizevector"  verbose = true begin
    eps_level = 1e-4

    @testset "Single vector" begin
        v = Vec{3,Float64}(0.0, 0.0, pi)
        n = normalizevector(v)
        @test n isa typeof(v)
        @test n == [0.0, 0.0, 1.0]

        v = Vec{3,Float64}(1.0, 1.0, 0.0)
        n = normalizevector(v)
        @test n isa typeof(v)
        @test isapprox(n,√2/2*[1.0, 1.0, 0.0], atol=eps_level)

        v = Point{3, Float64}( 1.0,  1.0,  1.0)
        n = normalizevector(v)
        @test n isa typeof(v)
        @test isapprox(n,√3/3*[1.0, 1.0, 1.0], atol=eps_level)
    end

    @testset "Vector set" begin
        M = cube(1.0)
        F = faces(M)
        V = coordinates(M)
        N = facenormal(F,V)
        U = [Vec{3,Float64}(n./(1.0.+rand(1))) for n in N]
        NN = normalizevector(U)

        @test eltype(NN) == eltype(U)
        @test isapprox(NN,N, atol=eps_level)
    end
end


@testset "circlepoints" verbose = true begin
    eps_level = 1e-4

    @testset "Errors" begin
        r = 2.0
        n = 40
        @test_throws ArgumentError circlepoints(2.0, n; dir=:wrong)
    end
    
    @testset "with value" begin
        r = 2.0
        n = 40
        V = circlepoints(r, n; dir=:acw)
        ind = round.(Int64,range(1,length(V),5))
        d = [sqrt(sum(v.^2)) for v in V]
        @test V isa Vector{Point3{Float64}}
        @test length(V) == n
        @test isapprox(d,fill(r,n),atol=eps_level)
        @test isapprox(V[ind], Point3{Float64}[[2.0, 0.0, 0.0], 
        [1.2246467991473532e-16, 2.0, 0.0], [-1.9753766811902753, 0.31286893008046196, 0.0], 
        [-0.31286893008046207, -1.9753766811902753, 0.0], 
        [1.9753766811902753, -0.31286893008046224, 0.0]], atol=eps_level)
   
        V = circlepoints(r, n; dir=:cw)
        ind = round.(Int64,range(1,length(V),5))
        d = [sqrt(sum(v.^2)) for v in V]
        @test V isa Vector{Point3{Float64}}
        @test length(V) == n
        @test isapprox(d,fill(r,n),atol=eps_level)
        @test isapprox(V[ind], Point3{Float64}[[2.0, 0.0, 0.0], [1.2246467991473532e-16, -2.0, 0.0], 
        [-1.9753766811902753, -0.31286893008046196, 0.0], [-0.31286893008046207, 1.9753766811902753, 0.0], 
        [1.9753766811902753, 0.31286893008046224, 0.0]], atol=eps_level)
    end

    @testset "with function" begin
        r = 1.0
        n = 40
        rFun(t) = r + 0.5 .* sin(3 * t)
        V = circlepoints(rFun, n)
        ind = round.(Int64,range(1,length(V),5))

        @test V isa Vector{Point3{Float64}}
        @test length(V) == n
        @test isapprox(V[ind],Point3{Float64}[[1.0, 0.0, 0.0], [3.061616997868383e-17, 0.5, 0.0],
        [-1.2118889022619925, 0.1919443455202825, 0.0], [-0.22612652951961248, -1.4277067182626626, 0.0], 
        [0.7634877789282826, -0.12092458456017953, 0.0]], atol=eps_level)
       
        V = circlepoints(rFun, n; dir=:cw)
        ind = round.(Int64,range(1,length(V),5))

        @test V isa Vector{Point3{Float64}}
        @test length(V) == n
        @test isapprox(V[ind],Point3{Float64}[[1.0, 0.0, 0.0], [9.184850993605148e-17, -1.5, 0.0], 
        [-0.7634877789282828, -0.12092458456017946, 0.0], [-0.08674240056084957, 0.5476699629276127, 0.0],
         [1.2118889022619928, 0.19194434552028272, 0.0]], atol=eps_level)
    end

end


@testset "loftlinear" verbose = true begin
    eps_level = 1e-4

    r = 1.0
    nc = 5
    n = Vec{3,Float64}(pi, pi, pi) # Offset vector
    V1 = circlepoints(r, nc; dir=:cw)
    V2 = [v.+n for v in V1]
    num_steps = 3

    @testset "Errors" begin
        @test_throws ArgumentError loftlinear(V1,V2;num_steps=num_steps,close_loop=true,face_type=:wrong)
    end

    @testset "quad" begin
        F,V = loftlinear(V1,V2;num_steps=num_steps,close_loop=true,face_type=:quad)
        ind = round.(Int64,range(1,length(V),5))

        @test F isa Vector{QuadFace{Int64}}
        @test length(F) == nc*(num_steps-1)
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test isapprox(V[ind], Point3{Float64}[[1.0, 0.0, 0.0], [-0.8090169943749475, 0.587785252292473, 0.0], [0.7617793324199491, 0.9830110745024232, 1.5707963267948966], [3.4506096479647406, 2.1905361372946395, 3.141592653589793], [3.4506096479647406, 4.092649169884947, 3.141592653589793]], atol=eps_level)

        # Providing "nothing" for the number of steps should create point spacing based number of steps
        F,V = loftlinear(V1,V2;num_steps=nothing,close_loop=true,face_type=:quad)
        @test length(V)/nc == ceil(Int64,norm(n)/(0.5*(pointspacingmean(V1)+pointspacingmean(V2))))

    end

    @testset "tri" begin
        F,V = loftlinear(V1,V2;num_steps=num_steps,close_loop=true,face_type=:tri)
        ind = round.(Int64,range(1,length(V),5))

        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == (nc*(num_steps-1))*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test isapprox(V[ind], Point3{Float64}[[1.0, 0.0, 0.0], [-0.8090169943749475, 0.587785252292473, 0.0], [0.7617793324199491, 1.5707963267948963, 1.5707963267948966], [3.4506096479647406, 2.1905361372946395, 3.141592653589793], [3.4506096479647406, 4.092649169884947, 3.141592653589793]], atol=eps_level)
    end

    @testset "tri_slash" begin
        F,V = loftlinear(V1,V2;num_steps=num_steps,close_loop=true,face_type=:tri_slash)
        ind = round.(Int64,range(1,length(V),5))

        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == (nc*(num_steps-1))*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test isapprox(V[ind], Point3{Float64}[[1.0, 0.0, 0.0], [-0.8090169943749475, 0.587785252292473, 0.0], [0.7617793324199491, 0.9830110745024232, 1.5707963267948966], [3.4506096479647406, 2.1905361372946395, 3.141592653589793], [3.4506096479647406, 4.092649169884947, 3.141592653589793]], atol=eps_level)
    end

    @testset "quad2tri" begin
        F,V = loftlinear(V1,V2;num_steps=num_steps,close_loop=true,face_type=:quad2tri)
        ind = round.(Int64,range(1,length(V),5))

        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == (nc*(num_steps-1))*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test isapprox(V[ind], Point3{Float64}[[1.0, 0.0, 0.0], [-0.8090169943749475, 0.587785252292473, 0.0], [0.7617793324199491, 0.9830110745024232, 1.5707963267948966], [3.4506096479647406, 2.1905361372946395, 3.141592653589793], [3.4506096479647406, 4.092649169884947, 3.141592653589793]], atol=eps_level)
    end
end

@testset "loftpoints2surf" verbose = true begin    

    nc = 5
    t = range(0,2.0*π-(2.0*π)/nc,nc) 
    V1 = [Point3{Float64}(cos(tt), sin(tt),0.0) for tt in t]
    V2 = [Point3{Float64}(v[1],v[2],v[3]+2.0) for v in V1]
    V3 = [Point3{Float64}(v[1],v[2],v[3]+4.0) for v in V1]
    Vp2 = vcat(V1,V2)
    Vp3 = vcat(V1,V2,V3)     

    @testset "errors" begin
        num_steps=4 # i.e. too many
        @test_throws ArgumentError  loftpoints2surf(Vp3,num_steps; close_loop=false,face_type=:quad)
    end    

    @testset "quad" begin
        #Non-closed loop
        num_steps = 2
        F,V = loftpoints2surf(Vp2,num_steps; close_loop=false,face_type=:quad)
        @test F isa Vector{QuadFace{Int64}}
        @test length(F) == (nc-1)*(num_steps-1)
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test V == Vp2

        # Closed loop
        num_steps = 2
        F,V = loftpoints2surf(Vp2,num_steps; close_loop=true,face_type=:quad)
        @test F isa Vector{QuadFace{Int64}}
        @test length(F) == nc*(num_steps-1)
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test V == Vp2

        # 3 layers
        num_steps = 3
        F,V = loftpoints2surf(Vp3,num_steps; close_loop=true,face_type=:quad)
        @test F isa Vector{QuadFace{Int64}}
        @test length(F) == nc*(num_steps-1)
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test V == Vp3
    end

    @testset "tri_slash" begin
        #Non-closed loop
        num_steps = 2
        F,V = loftpoints2surf(Vp2,num_steps; close_loop=false,face_type=:tri_slash)
        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == (nc-1)*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test V == Vp2

        # Closed loop
        num_steps = 2
        F,V = loftpoints2surf(Vp2,num_steps; close_loop=true,face_type=:tri_slash)
        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == nc*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test V == Vp2

        # 3 layers
        num_steps = 3
        F,V = loftpoints2surf(Vp3,num_steps; close_loop=true,face_type=:tri_slash)
        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == nc*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test V == Vp3
    end

    @testset "tri" begin
        #Non-closed loop
        num_steps = 2
        F,V = loftpoints2surf(Vp2,num_steps; close_loop=false,face_type=:tri)
        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == (nc-1)*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test V == Vp2

        # Closed loop
        num_steps = 2
        F,V = loftpoints2surf(Vp2,num_steps; close_loop=true,face_type=:tri)
        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == nc*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test V == Vp2

        # 3 layers
        num_steps = 3
        F,V = loftpoints2surf(Vp3,num_steps; close_loop=true,face_type=:tri)
        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == nc*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test V == Vp3
    end

    @testset "quad2tri" begin
        #Non-closed loop
        num_steps = 2
        F,V = loftpoints2surf(Vp2,num_steps; close_loop=false,face_type=:quad2tri)
        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == (nc-1)*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test V == Vp2

        # Closed loop
        num_steps = 2
        F,V = loftpoints2surf(Vp2,num_steps; close_loop=true,face_type=:quad2tri)
        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == nc*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test V == Vp2

        # 3 layers
        num_steps = 3
        F,V = loftpoints2surf(Vp3,num_steps; close_loop=true,face_type=:quad2tri)
        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == nc*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test V == Vp3
    end
end


@testset "dirplot" verbose = true begin
    M = cube(1.0)
    F = faces(M)
    V = coordinates(M)
    U = vertexnormal(F,V)

    fig = Figure(size=(800,800))
    ax = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Direction data plot")
    hp = poly!(ax,M, strokewidth=3,color=:white, shading = FastShading)
    hp1 = dirplot(ax,V,U; color=:black,linewidth=3,scaleval=1.0,style=:from)
    hp2 = dirplot(ax,V,U; color=:black,linewidth=3,scaleval=1.0,style=:to)
    hp3 = dirplot(ax,V,U; color=:black,linewidth=3,scaleval=1.0,style=:through)
    
    Mp = hp1[1].val

    @testset "Errors" begin
        @test_throws ArgumentError dirplot(ax,V,U; color=:black,linewidth=3,scaleval=1.0,style=:wrong)
    end

    @testset "Styles" begin
        @test typeof(hp1) == Wireframe{Tuple{GeometryBasics.Mesh{3, Float64, Line{3, Float64}, SimpleFaceView{3, Float64, 2, Int64, Point3{Float64}, LineFace{Int64}}}}}
        @test typeof(hp2) == Wireframe{Tuple{GeometryBasics.Mesh{3, Float64, Line{3, Float64}, SimpleFaceView{3, Float64, 2, Int64, Point3{Float64}, LineFace{Int64}}}}}
        @test typeof(hp3) == Wireframe{Tuple{GeometryBasics.Mesh{3, Float64, Line{3, Float64}, SimpleFaceView{3, Float64, 2, Int64, Point3{Float64}, LineFace{Int64}}}}}        
        @test length(faces(Mp)) == length(V)
    end
end


@testset "normalplot" verbose = true begin
    M = cube(1.0)
    F = faces(M)
    V = coordinates(M)

    fig = Figure(size=(800,800))
    ax = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Direction data plot")
    hp = poly!(ax,M, strokewidth=3,color=:white, shading = FastShading)


    @testset "Errors" begin
        @test_throws ArgumentError normalplot(ax,M; type_flag=:wrong, color=:black,linewidth=3,scaleval=nothing)
    end

    @testset "type_flag options" begin
        hp1 =  normalplot(ax,M; type_flag=:face, color=:black,linewidth=3,scaleval=nothing)
        Mp = hp1[1].val
        @test typeof(hp1) == Wireframe{Tuple{GeometryBasics.Mesh{3, Float64, Line{3, Float64}, SimpleFaceView{3, Float64, 2, Int64, Point3{Float64}, LineFace{Int64}}}}}
        @test length(faces(Mp)) == length(F)

        hp1 =  normalplot(ax,M; type_flag=:vertex, color=:black,linewidth=3,scaleval=nothing)
        Mp = hp1[1].val
        @test typeof(hp1) == Wireframe{Tuple{GeometryBasics.Mesh{3, Float64, Line{3, Float64}, SimpleFaceView{3, Float64, 2, Int64, Point3{Float64}, LineFace{Int64}}}}}
        @test length(faces(Mp)) == length(V)

        hp1 =  normalplot(ax,F,V; type_flag=:vertex, color=:black,linewidth=3,scaleval=nothing)
        Mp = hp1[1].val
        @test typeof(hp1) == Wireframe{Tuple{GeometryBasics.Mesh{3, Float64, Line{3, Float64}, SimpleFaceView{3, Float64, 2, Int64, Point3{Float64}, LineFace{Int64}}}}}
        @test length(faces(Mp)) == length(V)

        fileName_mesh = joinpath(comododir(),"assets","obj","spot_control_mesh.obj")
        Mn = load(fileName_mesh)   
        F = tofaces(faces(Mn))
        V = topoints(coordinates(Mn))

        hp1 =  normalplot(ax,F,V; type_flag=:vertex, color=:black,linewidth=3,scaleval=nothing)
        Mp = hp1[1].val
        @test typeof(hp1) == Wireframe{Tuple{GeometryBasics.Mesh{3, Float32, Line{3, Float32}, SimpleFaceView{3, Float32, 2, Int64, Point3{Float32}, LineFace{Int64}}}}}
        @test length(faces(Mp)) == length(V)

        hp1 =  normalplot(ax,Mn; type_flag=:face)
        Mp = hp1[1].val
        @test typeof(hp1) == Wireframe{Tuple{GeometryBasics.Mesh{3, Float32, Line{3, Float32}, SimpleFaceView{3, Float32, 2, Int64, Point3{Float32}, LineFace{Int64}}}}}
        @test length(faces(Mp)) == length(F)

        # Not supported yet
        # hp1 =  normalplot(ax,F,V; type_flag=:vertex, color=:black,linewidth=3,scaleval=nothing)
        # Mp = hp1[1].val
        # @test typeof(hp1) == Wireframe{Tuple{GeometryBasics.Mesh{3, Float64, Line{3, Float64}, SimpleFaceView{3, Float64, 2, Int64, Point3{Float64}, LineFace{Int64}}}}}
        # @test length(faces(Mp)) == length(V)

    end
end


@testset "wrapindex" verbose = true begin
    eps_level = 1e-4

    @testset "single value" begin
        n = 5
        m = 2
        @test wrapindex(1,n) == 1
        @test wrapindex(2,n) == 2
        @test wrapindex(n,n) == n
        @test wrapindex(n+m,n) == m      
    end

    @testset "Vector" begin
        n = 5
        a = 2
        b = 2*n
        m = [1,2,n,n+a,n+b]
        @test wrapindex(m,n) == [1,2,n,a,n]  
    end

    @testset "Unit range" begin
        n = 5
        a = 2
        b = 2*n
        m = 1:(2*n)+2
        r = wrapindex(m,n)
        @test r == [1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2]  
    end

    @testset "Step range" begin
        n = 5
        a = 2
        b = 2*n
        m = 1:2:(2*n)+2
        r = wrapindex(m,n)
        @test r == [1, 3, 5, 2, 4, 1] 
    end
end


@testset "edgeangles" begin
    eps_level = 1e-4
    # Regular cube
    M = cube(1.0)
    F = faces(M)
    V = coordinates(M)

    # Build deformation gradient tensor to induce shear with known angles
    fDef = zeros(3,3)
    for i=1:3
        fDef[i,i]=1.0
    end
    a = pi/4 # "45 degree shear"  
    fDef[1,2] = tan(a) 

    # Sheared cube coordinates
    V2 = topoints([fDef*v for v in V]) 

    A = edgeangles(F,V) # Angles for regular cube
    A2 = edgeangles(F,V2) # Angles for sheared cube
    @test all([all(a.==pi/2) for a in A]) # All right angles in undeformed cube
    @test isapprox(sort(unique(reduce(vcat,A2))),[pi/4, pi/2, pi/2+pi/4],atol=eps_level)
end


@testset "quad2tri" begin
    M = cube(1.0)
    F = faces(M)
    V = coordinates(M)

    # Build deformation gradient tensor to induce shear with known angles
    f = zeros(3,3)
    for i=1:3
        f[i,i]=1.0
    end
    a = pi/6 # "45 degree shear"  
    f[1,2] = tan(a) 
    V = topoints([f*v for v in V]) # Shear the cube

    # Build deformation gradient tensor to induce shear with known angles
    f = zeros(3,3)
    for i=1:3
        f[i,i]=1.0
    end
    a = pi/6 # "45 degree shear"  
    f[3,1] = tan(a) 
    V = topoints([f*v for v in V]) # Shear the cube

    @testset "Errors" begin
        @test_throws ArgumentError quad2tri(F,V; convert_method = :wrong)
    end

    @testset "convert_method variations" begin
        Ft = quad2tri(F,V; convert_method = :forward)
        @test Ft ==TriangleFace{Int64}[TriangleFace(1, 2, 3), TriangleFace(3, 4, 1), 
        TriangleFace(8, 7, 6), TriangleFace(6, 5, 8), TriangleFace(5, 6, 2), 
        TriangleFace(2, 1, 5), TriangleFace(6, 7, 3), TriangleFace(3, 2, 6), 
        TriangleFace(7, 8, 4), TriangleFace(4, 3, 7), TriangleFace(8, 5, 1), TriangleFace(1, 4, 8)]
        
        Ft = quad2tri(F,V; convert_method = :backward)
        @test Ft == TriangleFace{Int64}[TriangleFace(1, 2, 4), TriangleFace(2, 3, 4), 
        TriangleFace(8, 7, 5), TriangleFace(7, 6, 5), TriangleFace(5, 6, 1), 
        TriangleFace(6, 2, 1), TriangleFace(6, 7, 2), TriangleFace(7, 3, 2), 
        TriangleFace(7, 8, 3), TriangleFace(8, 4, 3), TriangleFace(8, 5, 4), TriangleFace(5, 1, 4)]

        Ft = quad2tri(F,V; convert_method = :angle)
        @test Ft == TriangleFace{Int64}[TriangleFace(1, 2, 4), TriangleFace(2, 3, 4), 
        TriangleFace(8, 7, 6), TriangleFace(6, 5, 8), TriangleFace(5, 6, 2), 
        TriangleFace(2, 1, 5), TriangleFace(6, 7, 3), TriangleFace(3, 2, 6), 
        TriangleFace(7, 8, 3), TriangleFace(8, 4, 3), TriangleFace(8, 5, 4), TriangleFace(5, 1, 4)]
    end
end


@testset "remove_unused_vertices" begin
    F, V = geosphere(2, 1.0)
    VC = simplexcenter(F, V)
    F = [F[i] for i in findall(map(v -> v[3] > 0, VC))] # Remove some faces
    Fn, Vn = remove_unused_vertices(F, V)

    @test Fn isa Vector{TriangleFace{Int64}}
    @test Vn isa Vector{Point3{Float64}}
    @test length(F) == length(Fn)
    @test length(Vn) ==  maximum(reduce(vcat,Fn))
    
    # Check empty results
    Fn, Vn = remove_unused_vertices(Vector{TriangleFace{Int64}}(), V)
    @test isempty(Fn)
    @test isempty(Vn)

end


@testset "trisurfslice" begin
    tol_level = 1e-2

    r = 2.5 # Sphere radius
    F,V = geosphere(3,r)

    p = [0.0,0.0,0.0]; # Point on cutting plane
    n = normalizevector(Vec{3, Float64}(0.0,0.0,1.0))# Cutting plane normal
    snapTolerance = 1e-6
    output_type = :full

    Fn,Vn,Cn = trisurfslice(F,V,n,p; output_type=output_type, snapTolerance=snapTolerance)

    # Check error
    @test_throws ArgumentError trisurfslice(F,V,n,p; output_type=:wrong)

    # Check if cut defines a circle of expected radius    
    Fn_below = Fn[Cn.<0]
    En_below = boundaryedges(Fn_below)
    ind_below = unique(reduce(vcat,En_below))
    d = [norm(v) for v in Vn[ind_below]]
    @test isapprox(sum((d.-r).^2),0.0,atol=tol_level) # Should 

    p = [0.0,0.0,0.0]; # Point on cutting plane
    n = normalizevector(Vec{3, Float64}(0.0,1.0,1.0))# Cutting plane normal
    snapTolerance = 1e-6
    output_type = :full

    Fn,Vn,Cn = trisurfslice(F,V,n,p; output_type=output_type)

    # Check if cut defines a circle of expected radius    
    Fn_below = Fn[Cn.<0]
    En_below = boundaryedges(Fn_below)
    ind_below = unique(reduce(vcat,En_below))
    d = [norm(v) for v in Vn[ind_below]]
    @test isapprox(sum((d.-r).^2),0.0,atol=tol_level) # Should 

    p = [0.0,0.0,0.0]; # Point on cutting plane
    n = normalizevector(Vec{3, Float64}(1.0,1.0,1.0))# Cutting plane normal
    snapTolerance = 1e-6
    output_type = :full

    Fn,Vn,Cn = trisurfslice(F,V,n,p; output_type=output_type)

    # Check if cut defines a circle of expected radius    
    Fn_below = Fn[Cn.<0]
    En_below = boundaryedges(Fn_below)
    ind_below = unique(reduce(vcat,En_below))
    d = [norm(v) for v in Vn[ind_below]]
    @test isapprox(sum((d.-r).^2),0.0,atol=tol_level) # Should 
end


@testset "count_edge_face" verbose = true begin

    
    @testset "Single triangle" begin
        F = [TriangleFace{Int64}(1, 2, 3)]       
        E = meshedges(F)
        Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        count_E2F = count_edge_face(F,Eu,indReverse)
        @test count_E2F == [1,1,1]

        # Reduced input set 
        count_E2F = count_edge_face(F)
        @test count_E2F == [1,1,1]

    end

    @testset "Single quad" begin
        F = [QuadFace{Int64}(1, 2, 3, 4)]       
        E = meshedges(F)
        Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        count_E2F = count_edge_face(F,Eu,indReverse)
        @test count_E2F == [1,1,1,1]
    end

    @testset "Triangles" begin
        F = [TriangleFace{Int64}(1, 2, 3),TriangleFace{Int64}(1, 4, 3)]
        E = meshedges(F)
        Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        count_E2F = count_edge_face(F,Eu,indReverse)
        @test count_E2F == [1,1,1,1,2]
    end

    @testset "Quads" begin
        F = [QuadFace{Int64}(1, 2, 3, 4),QuadFace{Int64}(6, 5, 4, 3)]
        E = meshedges(F)
        Eu,_,indReverse = gunique(E; return_unique=true, return_index=true, return_inverse=true, sort_entries=true)
        count_E2F = count_edge_face(F,Eu,indReverse)
        @test count_E2F == [1,1,1,1,2,1,1] 
    end

end


@testset "boundaryedges" verbose = true begin

    @testset "Single triangle" begin
        F = [TriangleFace{Int64}(1, 2, 3)]       
        Eb = boundaryedges(F)
        @test Eb == LineFace{Int64}[[1, 2], [2, 3], [3, 1]]
    end

    @testset "Single quad" begin
        F = [QuadFace{Int64}(1, 2, 3, 4)]       
        Eb = boundaryedges(F)
        @test Eb == LineFace{Int64}[[1, 2], [2, 3], [3, 4], [4, 1]]
    end

    @testset "Triangles" begin
        F = [TriangleFace{Int64}(1, 2, 3),TriangleFace{Int64}(1, 4, 3)]
        Eb = boundaryedges(F)
        @test Eb == LineFace{Int64}[[1, 2], [1, 4], [2, 3], [4, 3]]
    end

    @testset "Quads" begin
        F = [QuadFace{Int64}(1, 2, 3, 4),QuadFace{Int64}(6, 5, 4, 3)]
        Eb = boundaryedges(F)
        @test Eb == LineFace{Int64}[[1, 2], [6, 5], [2, 3], [5, 4], [4, 1], [3, 6]]
    end

end


@testset "edges2curve" begin

    # Open curve
    E = LineFace{Int64}[[1, 2], [2, 3], [3, 4], [4, 5]]
    ind = edges2curve(E)
    @test ind == [1,2,3,4,5]

    # Closed curve
    E = LineFace{Int64}[[1, 2], [2, 3], [3, 1]]
    ind = edges2curve(E)
    @test ind == [1,2,3,1]

end


@testset "pointspacingmean" verbose = true begin
    eps_level = 1e-4
    @testset "Curve" begin
        V = Point3{Float64}[[0.0,0.0,0.0],[0.25,0.0,0.0],[0.75,0.0,0.0],[1.75,0.0,0.0]]
        r = pointspacingmean(V)
        @test isapprox(r,mean(norm.(diff(V,dims=1))),atol = eps_level)
    end

    @testset "Edges" begin
        V = Point3{Float64}[[0.0,0.0,0.0],[0.25,0.0,0.0],[0.75,0.0,0.0],[1.75,0.0,0.0]]
        E = LineFace{Int64}[[1,2],[2,3],[3,4]]
        r = pointspacingmean(E,V)
        @test isapprox(r,mean(norm.(diff(V,dims=1))),atol = eps_level)
    end

    @testset "Faces" begin
        V = Point3{Float64}[[0.0,0.0,0.0],[0.25,0.0,0.0],[0.25,0.5,0.0],[0,0.5,0.0],[0.0,0.0,0.0]]
        F = QuadFace{Int64}[[1,2,3,4]]
        r = pointspacingmean(F,V)
        @test isapprox(r,mean(norm.(diff(V,dims=1))),atol = eps_level)
    end

    @testset "Mesh" begin        
        M = cube(sqrt(3))
        F = faces(M)
        V = coordinates(M)
        @test pointspacingmean(F,V)==2.0
        @test pointspacingmean(M)==2.0
    end
end


@testset "extrudecurve" verbose = true begin
    eps_level = 1e-4

    r = 1.0
    nc = 16
    d = 3.0
    Vc = circlepoints(r, nc; dir=:cw)
    num_steps = 5

    @testset "Default behaviours" begin
        F, V = extrudecurve(Vc, d)
        z = [v[3] for v in V]
        zMax = maximum(z)
        zMin = minimum(z)            
        @test isapprox(zMax,d,atol = eps_level) && isapprox(zMin,0.0,atol = eps_level)

        F, V = extrudecurve(Vc, d; face_type=:tri)
        z = [v[3] for v in V]
        zMax = maximum(z)
        zMin = minimum(z)            
        @test isapprox(zMax,d,atol = eps_level) && isapprox(zMin,0.0,atol = eps_level)
    end

    @testset "Direction (s) variations" begin
        F, V = extrudecurve(Vc, d; s=1, n=Vec{3, Float64}(0.0,0.0,1.0), num_steps=num_steps, close_loop=true, face_type=:quad)
        z = [v[3] for v in V]
        zMax = maximum(z)
        zMin = minimum(z)            
        @test isapprox(zMax,d,atol = eps_level) && isapprox(zMin,0.0,atol = eps_level)
       
        F, V = extrudecurve(Vc, d; s=0, n=Vec{3, Float64}(0.0,0.0,1.0), num_steps=num_steps, close_loop=true, face_type=:quad)
        z = [v[3] for v in V]
        zMax = maximum(z)
        zMin = minimum(z)            
        @test isapprox(zMax,d/2,atol = eps_level) && isapprox(zMin,-d/2,atol = eps_level)

        F, V = extrudecurve(Vc, d; s=-1, n=Vec{3, Float64}(0.0,0.0,1.0), num_steps=num_steps, close_loop=true, face_type=:quad)
        z = [v[3] for v in V]
        zMax = maximum(z)
        zMin = minimum(z)            
        @test isapprox(zMax,0.0,atol = eps_level) && isapprox(zMin,-d,atol = eps_level)
    end
    
    @testset "Direction (n) variations" begin
        n=Vec{3, Float64}(0.0,0.0,1.0) # Upward
        F, V = extrudecurve(Vc, d; s=1, n=n, num_steps=num_steps, close_loop=true, face_type=:quad)
        z = [v[3] for v in V]
        zMax = maximum(z)
        zMin = minimum(z)            
        @test isapprox(zMax,d,atol = eps_level) && isapprox(zMin,0.0,atol = eps_level)
       
        n=Vec{3, Float64}(0.0,0.0,-1.0) # Downward
        F, V = extrudecurve(Vc, d; s=1, n=n, num_steps=num_steps, close_loop=true, face_type=:quad)
        z = [v[3] for v in V]
        zMax = maximum(z)
        zMin = minimum(z)            
        @test isapprox(zMax,0.0,atol = eps_level) && isapprox(zMin,-d,atol = eps_level)

        n = normalizevector(Vec{3, Float64}(1.0,0.0,1.0)) # 45 degree direction upward
        F, V = extrudecurve(Vc, d; s=1, n=n, num_steps=num_steps, close_loop=true, face_type=:quad)
        z = [v[3] for v in V]
        zMax = maximum(z)
        zMin = minimum(z)            
        @test isapprox(zMax,sqrt(2.0)*d/2,atol = eps_level) && isapprox(zMin,0.0,atol = eps_level)
    end

    @testset "face_type=:quad" begin
        F, V = extrudecurve(Vc, d; s=1, n=Vec{3, Float64}(0.0,0.0,1.0), num_steps=num_steps, close_loop=true, face_type=:quad)
        z = [v[3] for v in V]
        zMax = maximum(z)
        zMin = minimum(z)

        @test F isa Vector{QuadFace{Int64}}
        @test length(F) == nc*(num_steps-1)

        ind = round.(Int64,range(1,length(V),5))
        @test V isa Vector{Point3{Float64}}
        @test isapprox(zMax,d,atol = eps_level) && isapprox(zMin,0.0,atol = eps_level)
        @test isapprox(V[ind],Point3{Float64}[[1.0, 0.0, 0.0], [2.83276944882399e-16, -1.0, 0.75], 
        [-0.9238795325112867, -0.3826834323650899, 1.5], [-0.38268343236509034, 0.9238795325112865, 2.25], 
        [0.9238795325112865, 0.3826834323650904, 3.0]],atol = eps_level)
    end

    @testset "face_type=:tri" begin
        F, V = extrudecurve(Vc, d; s=1, n=Vec{3, Float64}(0.0,0.0,1.0), num_steps=num_steps, close_loop=true, face_type=:tri)
        z = [v[3] for v in V]
        zMax = maximum(z)
        zMin = minimum(z)

        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == (nc*(num_steps-1))*2

        ind = round.(Int64,range(1,length(V),5))
        @test V isa Vector{Point3{Float64}}
        @test isapprox(zMax,d,atol = eps_level) && isapprox(zMin,0.0,atol = eps_level)
        @test isapprox(V[ind],Point3{Float64}[[1.0, 0.0, 0.0], [-0.1913417161825446, -0.9619397662556435, 0.75], 
        [-0.9238795325112867, -0.3826834323650899, 1.5], [-0.19134171618254525, 0.9619397662556433, 2.25], 
        [0.9238795325112865, 0.3826834323650904, 3.0]],atol = eps_level)
    end

    @testset "face_type=:tri_slash" begin
        F, V = extrudecurve(Vc, d; s=1, n=Vec{3, Float64}(0.0,0.0,1.0), num_steps=num_steps, close_loop=true, face_type=:tri_slash)
        z = [v[3] for v in V]
        zMax = maximum(z)
        zMin = minimum(z)

        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == (nc*(num_steps-1))*2

        ind = round.(Int64,range(1,length(V),5))
        @test V isa Vector{Point3{Float64}}
        @test isapprox(zMax,d,atol = eps_level) && isapprox(zMin,0.0,atol = eps_level)
        @test isapprox(V[ind],Point3{Float64}[[1.0, 0.0, 0.0], [2.83276944882399e-16, -1.0, 0.75], 
        [-0.9238795325112867, -0.3826834323650899, 1.5], [-0.38268343236509034, 0.9238795325112865, 2.25], 
        [0.9238795325112865, 0.3826834323650904, 3.0]],atol = eps_level)
    end

    @testset "face_type=:quad2tri" begin
        F, V = extrudecurve(Vc, d; s=1, n=Vec{3, Float64}(0.0,0.0,1.0), num_steps=num_steps, close_loop=true, face_type=:quad2tri)
        z = [v[3] for v in V]
        zMax = maximum(z)
        zMin = minimum(z)

        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == (nc*(num_steps-1))*2

        ind = round.(Int64,range(1,length(V),5))
        @test V isa Vector{Point3{Float64}}
        @test isapprox(zMax,d,atol = eps_level) && isapprox(zMin,0.0,atol = eps_level)
        @test isapprox(V[ind],Point3{Float64}[[1.0, 0.0, 0.0], [2.83276944882399e-16, -1.0, 0.75], 
        [-0.9238795325112867, -0.3826834323650899, 1.5], [-0.38268343236509034, 0.9238795325112865, 2.25], 
        [0.9238795325112865, 0.3826834323650904, 3.0]],atol = eps_level)
    end

end


@testset "meshgroup" verbose = true begin

    @testset "Errors" begin
        F = TriangleFace{Int64}[[1,2,3]]
        @test_throws ArgumentError meshgroup(F; con_type=:wrong)
    end

    @testset "Single face" begin
        # Single triangle
        F = TriangleFace{Int64}[[1,2,3]]
        C = meshgroup(F)
        @test C == [1]

        # Single triangle, edge connectivity
        F = TriangleFace{Int64}[[1,2,3]]
        C = meshgroup(F; con_type=:e)
        @test C == [1]
        
        # Single quad
        F = QuadFace{Int64}[[1,2,3,4]]
        C = meshgroup(F)
        @test C == [1]
    end

    @testset "Two separate faces" begin
        # Two triangles separate
        F = TriangleFace{Int64}[[1,2,3],[4,5,6]]
        C = meshgroup(F)
        @test C == [1,2]

         # Two triangles, edge sharing 
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]
        C = meshgroup(F; con_type=:e)
        @test C == [1,1]
        
        F = TriangleFace{Int64}[[1,2,3],[2,3,4]]
        C = meshgroup(F; con_type=:v)
        @test C == [1,1]

        # Two bow-tie triangles, node sharing 
        F = TriangleFace{Int64}[[1,2,3],[3,4,5]]
        C = meshgroup(F; con_type=:e)
        @test C == [1,2]

        F = TriangleFace{Int64}[[1,2,3],[3,4,5]]
        C = meshgroup(F; con_type=:v)
        @test C == [1,1]
        
        # Two quads
        F = QuadFace{Int64}[[1,2,3,4],[5,6,7,8]]
        C = meshgroup(F)
        @test C == [1,2]
    end

    @testset "Single group" begin
        # Single tetrahedron
        M = tetrahedron(1.0)
        F = faces(M)
        C = meshgroup(F)
        @test C == ones(length(F))

        # Single cube
        M = cube(1.0)
        F = faces(M)
        C = meshgroup(F)
        @test C == ones(length(F))
    end

    @testset "Two groups" begin
        # Three triangles two joined, one separate
        F = TriangleFace{Int64}[[1,2,3],[2,3,4],[5,6,7]]
        C = meshgroup(F)
        @test C == [1,1,2]
        
        # Two tetrahedrons
        M = tetrahedron(1.0)
        F = faces(M)
        V = coordinates(M)
        n = length(F)
        F2 = map(f-> f.+length(V),F)
        V2 = map(v-> Point{3, Float64}(2.0+v[1],v[2],v[3]),V)
        append!(F,F2)
        append!(V,V2)
        C = meshgroup(F)
        @test C == repeat(1:2,inner=n)

        # Two tetrahedrons
        M = cube(1.0)
        F = faces(M)
        V = coordinates(M)
        n = length(F)
        F2 = map(f-> f.+length(V),F)
        V2 = map(v-> Point{3, Float64}(2.0+v[1],v[2],v[3]),V)
        append!(F,F2)
        append!(V,V2)
        C = meshgroup(F)
        @test C == repeat(1:2,inner=n)
    end

    
end


@testset "distmarch" verbose=true begin
    eps_level = 1e-4

    @testset "Single face" begin
        # Single triangle
        F = TriangleFace{Int64}[[1,2,3]]
        V = Point3{Float64}[[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0]]
        d,dd,l = distmarch(F,V,[1])
        @test isapprox(d,[0.0,1.0,sqrt(2)],atol=eps_level)

        # Single triangle, un-used nodes
        F = TriangleFace{Int64}[[1,2,3]]
        V = Point3{Float64}[[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0],
                            [1.0,0.0,0.0],[1.0,1.0,0.0]]
        d,dd,l = distmarch(F,V,[1])
        r = [0.0,1.0,sqrt(2),NaN,NaN]
        b = .!isnan.(r)
        @test isapprox(d[b],r[b],atol=eps_level) # Check reached entries
        @test all(isnan.(d[.!b])) # Now check NaNs

        # Single quad
        F = QuadFace{Int64}[[1,2,3,4]]
        V = Point3{Float64}[[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0],[0.0,1.0,0.0]]
        d,dd,l = distmarch(F,V,[1])
        @test isapprox(d,[0.0,1.0,sqrt(2),1.0],atol=eps_level)
    end

    @testset "Multi-face meshes" begin
        # Bowtie triangle set
        V = Point3{Float64}[[0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0], 
                            [2.0, 1.0, 0.0], [2.0, -1.0, 0.0]]
        F = TriangleFace{Int64}[TriangleFace(1, 2, 3), TriangleFace(5, 4, 3)]
        d,dd,l = distmarch(F,V,[1])
        @test isapprox(d,[0.0,2.0,sqrt(2),2*sqrt(2),2*sqrt(2)],atol=eps_level)
        
        # Two disconnected triangles
        V = Point3{Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], 
                            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]]
        F = TriangleFace{Int64}[TriangleFace(1, 2, 3), TriangleFace(4, 5, 6)]
        d,dd,l = distmarch(F,V,[1])
        r = [0.0,1.0,sqrt(2),NaN,NaN,NaN]
        b = .!isnan.(r)
        @test isapprox(d[b],r[b],atol=eps_level) # Check reached entries
        @test all(isnan.(d[.!b])) # Now check NaNs

        # Two disconnected quads
        V = Point3{Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0,1.0,0.0],
                            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0,1.0,1.0]]
        F = QuadFace{Int64}[QuadFace(1, 2, 3, 4), QuadFace(5, 6, 7, 8)]
        d,dd,l = distmarch(F,V,[1])
        r = [0.0,1.0,sqrt(2),1.0,NaN,NaN,NaN,NaN]
        b = .!isnan.(r)
        @test isapprox(d[b],r[b],atol=eps_level) # Check reached entries
        @test all(isnan.(d[.!b])) # Now check NaNs

        # Single cube
        r = sqrt(3)
        M = cube(r)
        F = faces(M)
        V = coordinates(M)
        d,dd,l = distmarch(F,V,[1])
        @test isapprox(d,[0.0, 2.0, 2.8284271247461903, 2.0, 2.0,
                         2.8284271247461903, 4.82842712474619, 2.8284271247461903],atol=eps_level) 
        
        # Triangulated sphere, distance should approximate π 
        r = 1.0
        F,V = geosphere(4,r)
        z = [v[3] for v in V]
        indStart =[findmin(z)[2]]
        d,dd,l = distmarch(F,V,indStart)
        @test isapprox(maximum(d),π,atol=0.01)

        # Quadrangulated sphere, distance should approximate π 
        r = 1.0
        F,V = quadsphere(4,r)
        z = [v[3] for v in V]
        indStart =[findmin(z)[2]]
        d,dd,l = distmarch(F,V,indStart)
        @test isapprox(maximum(d),π,atol=0.01)
    end
end


# @testset "distseedpoints" verbose = true begin

# end


@testset "ray_triangle_intersect" verbose = true begin
    eps_level = 1e-4

    # Single cube
    r = sqrt(3)
    M = cube(r)
    F = faces(M)
    V = coordinates(M)
    F = quad2tri(F,V,convert_method = :forward)
    @testset "ray" begin 
        ray_origin = GeometryBasics.Point3{Float64}(0.25,0.0,1.5) # Slight off so we hit one triangle, not two at the edge
        ray_vector = Vec3{Float64}(0.0,0.0,-1)

        P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = 1)
        @test isapprox(P,Point3{Float64}[[0.25, 0.0, 1.0]],atol=eps_level)
        @test isa(indIntersect,Vector{Int64}) # indIntersect == 3

        P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = 0)
        @test isapprox(P,Point3{Float64}[[0.25, 0.0, -1.0], [0.25, 0.0, 1.0]],atol=eps_level)
        @test isa(indIntersect,Vector{Int64})

        P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = -1)
        @test isapprox(P,Point3{Float64}[[0.25, 0.0, -1.0]],atol=eps_level)
        @test isa(indIntersect,Vector{Int64})

        ray_origin = GeometryBasics.Point3{Float64}(0.0,0.0,1.5) # At centre so hits an edge between two triangles
        ray_vector = Vec3{Float64}(0.0,0.0,-1)
        P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :ray, triSide = 0, tolEps = 1e-3)
        @test isapprox(P,Point3{Float64}[[0.0, 0.0, -1.0], [0.0, 0.0, -1.0], 
                                         [0.0, 0.0,  1.0], [0.0, 0.0,  1.0]],atol=eps_level)
        @test isa(indIntersect,Vector{Int64})
    end

    @testset "line type" begin 
        ray_origin = GeometryBasics.Point3{Float64}(0.25,0.0,1.5) # Slight off so we hit one triangle, not two at the edge
        ray_vector = Vec3{Float64}(0.0,0.0,-1) # Shorst so only one hit

        P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = 1)
        @test isapprox(P,Point3{Float64}[[0.25, 0.0, 1.0]],atol=eps_level)
        @test isa(indIntersect,Vector{Int64}) # indIntersect == 3

        ray_vector = Vec3{Float64}(0.0,0.0,-3) # Long so two hits potentially

        P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = 1)
        @test isapprox(P,Point3{Float64}[[0.25, 0.0, 1.0]],atol=eps_level)
        @test isa(indIntersect,Vector{Int64}) # indIntersect == 3

        P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = 0)
        @test isapprox(P,Point3{Float64}[[0.25, 0.0, -1.0], [0.25, 0.0, 1.0]],atol=eps_level)
        @test isa(indIntersect,Vector{Int64})

        P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = -1)
        @test isapprox(P,Point3{Float64}[[0.25, 0.0, -1.0]],atol=eps_level)
        @test isa(indIntersect,Vector{Int64})

        ray_origin = GeometryBasics.Point3{Float64}(0.0,0.0,1.5) # At centre so hits an edge between two triangles
        P,indIntersect = ray_triangle_intersect(F,V,ray_origin,ray_vector; rayType = :line, triSide = 0, tolEps = 1e-3)
        @test isapprox(P,Point3{Float64}[[0.0, 0.0, -1.0], [0.0, 0.0, -1.0], 
                                         [0.0, 0.0,  1.0], [0.0, 0.0,  1.0]],atol=eps_level)
        @test isa(indIntersect,Vector{Int64})
    end
end


@testset "mesh_curvature_polynomial" verbose = true begin
    tol_level = 0.01

    # A sphere has constant curvature. Both K1 and K2 are equivalent. curvature is 1/r 
    @testset "Triangulated sphere" begin
        r = 10 
        k_true = 1.0/r
        F,V = geosphere(3,r) 
        K1,K2,U1,U2,H,G = mesh_curvature_polynomial(F,V)
        @test isapprox(mean(K1),k_true,atol=tol_level)
        @test isapprox(mean(K2),k_true,atol=tol_level)
        @test isapprox(mean(H),k_true,atol=tol_level)
        @test isapprox(sqrt(abs(mean(G))),k_true,atol=tol_level)
    end

    @testset "Quadrangulated sphere" begin
        r = 10
        k_true = 1.0/r
        F,V = quadsphere(4,r) 
        K1,K2,U1,U2,H,G = mesh_curvature_polynomial(F,V)
        @test isapprox(mean(K1),k_true,atol=tol_level)
        @test isapprox(mean(K2),k_true,atol=tol_level)
        @test isapprox(mean(H),k_true,atol=tol_level)
        @test isapprox(sqrt(abs(mean(G))),k_true,atol=tol_level)
    end

    # For this algorithm the cube appears to have homogeneous curvature too that is 1.5/r
    @testset "Cube" begin
        r = sqrt(3) 
        k_true = 1.5*(1.0/r) # Note only for the polynomial method
        M = cube(r)
        F = faces(M)
        V = coordinates(M)
        K1,K2,U1,U2,H,G = mesh_curvature_polynomial(F,V)
        @test isapprox(mean(K1),k_true,atol=tol_level)
        @test isapprox(mean(K2),k_true,atol=tol_level)
        @test isapprox(mean(H),k_true,atol=tol_level)
        @test isapprox(sqrt(abs(mean(G))),k_true,atol=tol_level)

        # Check if mesh input functions the same
        K1m,K2m,U1m,U2m,Hm,Gm = mesh_curvature_polynomial(M)
        @test K1==K1m
        @test K2==K2m       
    end

    # For a cylinder k1=1/r while k2=0, hence H = 1/(2*R)
    @testset "Cylinder" begin
        r = 25
        s = r/10
        nc = round(Int64,(2*pi*r)/s)
        d = r*2
        Vc = circlepoints(r, nc; dir=:cw)
        num_steps = round(Int64,d/s)
        num_steps = num_steps + Int64(iseven(num_steps))
        F, V = extrudecurve(Vc, d; s=1, n=[0.0,0.0,1.0], num_steps=num_steps, close_loop=true, face_type=:quad)
        K1,K2,U1,U2,H,G = mesh_curvature_polynomial(F,V)
        @test isapprox(mean(K1),1.0/r,atol=tol_level)
        @test isapprox(mean(K2),0.0,atol=tol_level)
        @test isapprox(mean(H),1.0/(2*r),atol=tol_level)
        @test isapprox(sqrt(abs(mean(G))),0.0,atol=tol_level)
    end
end



@testset "separate_vertices" verbose = true begin
    
    @testset "Single face" begin
        # Single triangle
        F = TriangleFace{Int64}[[1,2,3]]
        V = Point3{Float64}[[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0]]
        Fn, Vn = separate_vertices(F, V)    
        @test Vn isa Vector{Point3{Float64}}
        @test typeof(Fn) == typeof(F)
        @test length(Fn) == length(F)    
        @test length(Vn) == length(F)*length(F[1])

        # Single triangle, un-used nodes
        F = TriangleFace{Int64}[[1,2,3]]
        V = Point3{Float64}[[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0],
                            [1.0,0.0,0.0],[1.0,1.0,0.0]]
        Fn, Vn = separate_vertices(F, V)    
        @test Vn isa Vector{Point3{Float64}}
        @test typeof(Fn) == typeof(F)
        @test length(Fn) == length(F)    
        @test length(Vn) == length(F)*length(F[1])

        # Single quad
        F = QuadFace{Int64}[[1,2,3,4]]
        V = Point3{Float64}[[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0],[0.0,1.0,0.0]]
        Fn, Vn = separate_vertices(F, V)    
        @test Vn isa Vector{Point3{Float64}}
        @test typeof(Fn) == typeof(F)
        @test length(Fn) == length(F)    
        @test length(Vn) == length(F)*length(F[1])
    end

    @testset "Quadrilateral face mesh" begin
        M = cube(1.0)
        F = faces(M)
        V = coordinates(M)
        Fn, Vn = separate_vertices(F, V)    
        @test Vn isa Vector{Point3{Float64}}
        @test typeof(Fn) == typeof(F)
        @test length(Fn) == length(F)    
        @test length(Vn) == length(F)*length(F[1])
    end

    @testset "Triangulated mesh" begin
        M = tetrahedron(1.0)
        F = faces(M)
        V = coordinates(M)
        Fn, Vn = separate_vertices(F, V)    
        @test Vn isa Vector{Point3{Float64}}
        @test typeof(Fn) == typeof(F)
        @test length(Fn) == length(F)    
        @test length(Vn) == length(F)*length(F[1])
    end

    @testset "Mesh" begin
        M = tetrahedron(1.0)
        F = faces(M)
        V = coordinates(M)
        Mn = separate_vertices(M)    
        Fn = faces(Mn)
        Vn = coordinates(Mn)
        @test Vn isa Vector{Point3{Float64}}
        @test typeof(Fn) == typeof(F)
        @test length(Fn) == length(F)    
        @test length(Vn) == length(F)*length(F[1])
    end

end


@testset "curve_length" verbose = true begin
    tol_level = 1e-6
    r = 2.25
    nc = 10    
    V = circlepoints(r, nc; dir=:cw)    
    length_true  = collect(range(0,(nc*2.0*r*sin(0.5*((2.0*pi)/nc))),nc+1))    
    L = curve_length(V; close_loop=false)
    @test isapprox(L,length_true[1:end-1],atol=tol_level)    

    L = curve_length(V; close_loop=true)
    @test isapprox(L,length_true,atol=tol_level)    
    @test L isa Vector{Float64}
end


@testset "evenly_sample" begin
    
    tol_level = 1e-4

    # Even sampling should be nearly perfect for a linear curve 
    @testset "Evenly upsampling linear curve" begin
        # Example linear curve raw
        X = range(0, 10, 5)
        V = [GeometryBasics.Point{3,Float64}(x, 2.0*x, 0.0) for x in X]

        # Create true evenly upsampled data 
        n = 20
        Xt = range(0, 10, n)
        Vt = [GeometryBasics.Point{3,Float64}(x, 2.0*x, 0.0) for x in Xt]

        # Resample using evenly_sample
        Vi, S = evenly_sample(V, n; niter=5)

        @test sum(norm.(Vi-Vt)) < tol_level # Correct and even spacing
        @test typeof(V) == typeof(Vi) # Did not manipulate input type
        @test length(Vi) == n # Correct length
    end

    # Even sampling should be nearly perfect for downsampling
    @testset "Evenly downsampling circular curve" begin
        # Example circle curve raw
        r = 2.25
        nc = 100    
        V = [GeometryBasics.Point{3, Float64}(r*cos(t),r*sin(t),0) for t in range(0.0,2.0*π,nc)]

        # Create true evenly upsampled data 
        n = 10
        Vt = [GeometryBasics.Point{3, Float64}(r*cos(t),r*sin(t),0) for t in range(0.0,2.0*π,n)]
        
        # Resample using evenly_sample
        Vi, S = evenly_sample(V, n; niter=5)

        @test sum(norm.(Vi-Vt)) < tol_level # Correct and even spacing
        @test typeof(V) == typeof(Vi) # Did not manipulate input type
        @test length(Vi) == n # Correct length
    end
end


@testset "evenly_space" verbose = true begin
    eps_level = 1e-3
    np = 10
    t = range(0.0,2.0*π,np) # Parameterisation metric
    V = [GeometryBasics.Point{3, Float64}(cos(t[i]),sin(t[i]),0.0) for i in eachindex(t)] 
    
    Vn = evenly_space(V)
    @test isapprox(pointspacingmean(Vn),pointspacingmean(V),atol=eps_level)

    pointSpacing = pointspacingmean(V)
    Vn = evenly_space(V,pointSpacing)
    @test isapprox(pointspacingmean(Vn),pointSpacing,atol=eps_level)
end


@testset "invert_faces" begin
    # Single face
    F = TriangleFace{Int64}[[1,2,3]]
    F_inv = TriangleFace{Int64}[[3,2,1]]
    @test invert_faces(F)==F_inv

    # Two face
    F = [TriangleFace{Int64}(1, 2, 3),TriangleFace{Int64}(4, 5, 6)]   
    F_inv = [TriangleFace{Int64}(3, 2, 1),TriangleFace{Int64}(6, 5, 4)]   
    @test invert_faces(F)==F_inv    
end

@testset "kabsch_rot" begin
    tol_level = 1e-6

    M = cube(sqrt(3))
    V1 = coordinates(M)
    R_true = RotXYZ(0.25*π,0.25*π,0.25*π)
    V2 = [R_true*v for v in V1]
    R_kabsch_forward = kabsch_rot(V1,V2)
    R_kabsch_backward = kabsch_rot(V2,V1)
    V1r = [R_kabsch_backward*v for v in V2]
    @test isapprox(R_kabsch_forward,R_true,atol=tol_level) # Able to retrieve forward rotation
    @test isapprox(V1r,V1,atol=tol_level) # Check if backward rotation is succesful   
end

@testset "sweeploft" verbose = true begin
    tol_level = 1e-6

    # Define guide curve
    nc = 25 # Number of points on guide curve
    P = Vector{GeometryBasics.Point{3, Float64}}(undef,4)
    P[1 ] = GeometryBasics.Point{3, Float64}( 0.0, 0.0, 0.0)
    P[2 ] = GeometryBasics.Point{3, Float64}( 1.0, 0.0, 0.0)
    P[3 ] = GeometryBasics.Point{3, Float64}( 1.0, 1.0, 0.0)
    P[4 ] = GeometryBasics.Point{3, Float64}( 1.0, 1.0, 1.0)
    Vc = nbezier(P,nc) # Get Bezier fit points
    Vc = [vc.*10 for vc in Vc]
    Vc,Sc = evenly_sample(Vc, nc)

    # Define section curves
    np = 20 # Number of section points    
    V1 = circlepoints(2.0,np; dir=:acw)
    V1,_ = evenly_sample(V1, np)
    Q = RotXYZ(0.0,0.5*π,0.0) # Define a rotation tensor using Euler angles
    V1 = [Q*v for v in V1] # Rotate the coordinates

    V2 = circlepoints(2.0,np; dir=:acw)
    V2,_ = evenly_sample(V2, np)
    V2 = [v2 .+ Vc[end] for v2 in V2] 

    @testset "reversed" begin
        F,V = sweeploft(Vc,V1,V2; face_type=:quad, num_twist=0, close_loop=true)
        Fr,Vr = sweeploft(Vc,reverse(V1),reverse(V2); face_type=:quad, num_twist=0, close_loop=true)

        @test F == Fr
        @test isa(F,Vector{QuadFace{Int64}})
        @test length(Vr) == nc*np
        @test length(Vr) == length(V)
        @test isa(V,Vector{Point3{Float64}})  
    end

    @testset "quad" begin
        F,V = sweeploft(Vc,V1,V2; face_type=:quad, num_twist=0, close_loop=true)
        @test length(F) == (nc-1)*np
        @test isa(F,Vector{QuadFace{Int64}})
        @test length(V) == nc*np
        @test isa(V,Vector{Point3{Float64}})

        F,V = sweeploft(Vc,V1,V2; face_type=:quad, num_twist=1, close_loop=false)
        @test length(F) == (nc-1)*(np-1)
        @test isa(F,Vector{QuadFace{Int64}})
        @test length(V) == nc*np
        @test isa(V,Vector{Point3{Float64}})        
    end

    @testset "tri_slash" begin
        F,V = sweeploft(Vc,V1,V2; face_type=:tri_slash, num_twist=0, close_loop=true)
        @test length(F) == (nc-1)*np*2
        @test isa(F,Vector{TriangleFace{Int64}})
        @test length(V) == nc*np
        @test isa(V,Vector{Point3{Float64}})

        F,V = sweeploft(Vc,V1,V2; face_type=:tri_slash, num_twist=1, close_loop=false)
        @test length(F) == (nc-1)*(np-1)*2
        @test isa(F,Vector{TriangleFace{Int64}})
        @test length(V) == nc*np
        @test isa(V,Vector{Point3{Float64}})        
    end

    @testset "tri" begin
        F,V = sweeploft(Vc,V1,V2; face_type=:tri, num_twist=0, close_loop=true)
        @test length(F) == (nc-1)*np*2
        @test isa(F,Vector{TriangleFace{Int64}})
        @test length(V) == nc*np
        @test isa(V,Vector{Point3{Float64}})

        F,V = sweeploft(Vc,V1,V2; face_type=:tri, num_twist=1, close_loop=false)
        @test length(F) == (nc-1)*(np-1)*2
        @test isa(F,Vector{TriangleFace{Int64}})
        @test length(V) == nc*np
        @test isa(V,Vector{Point3{Float64}})        
    end

    @testset "quad2tri" begin
        F,V = sweeploft(Vc,V1,V2; face_type=:quad2tri, num_twist=0, close_loop=true)
        @test length(F) == (nc-1)*np*2
        @test isa(F,Vector{TriangleFace{Int64}})
        @test length(V) == nc*np
        @test isa(V,Vector{Point3{Float64}})

        F,V = sweeploft(Vc,V1,V2; face_type=:quad2tri, num_twist=1, close_loop=false)
        @test length(F) == (nc-1)*(np-1)*2
        @test isa(F,Vector{TriangleFace{Int64}})
        @test length(V) == nc*np
        @test isa(V,Vector{Point3{Float64}})        
    end
end


@testset "revolvecurve" verbose = true begin
    tol_level = 1e-6
    nc = 5
    t = range(1.0,2.0,nc)
    Vc = [Point3(tt,0.0,0.0) for tt in t]
    n = Vec{3, Float64}(0.0,0.0,1.0)
    rMax = 0.0
    rMin = Inf
    for v in Vc
        typeof(n)
        typeof(v)
        rNow = dot(normalizevector(cross(cross(n,v),n)),v)
        if !isnan(rNow)
            rMax = max(rMax,rNow)
            rMin = min(rMin,rNow)
        else
            rMin = 0.0
        end
    end    
        
    num_steps = 6
    close_loop = false

    ind = round.(Int64,range(1,nc*num_steps,5))

    θ = 1.0*π

    # Check if revolved line which is orthogonal to n yields zero dot product 
    # with n, and check if this is rotation invariant. 
    @testset "Rotation invariance" begin
        _,VN = geosphere(3,1.0)   # Create point on a sphere
        b = zeros(Bool,length(VN))
        for i in eachindex(VN)
            m = VN[i] # Current rotation axis
            Q = rotation_between(n,m) # Rotation matrix to rotate Vc orthogonal to m       
            F1,V1 = revolvecurve([Q*v for v in Vc],θ;  s=1, n=m,num_steps=num_steps,close_loop=close_loop,face_type=:quad)
            
            d = sum([dot(m,v) for v in V1])
            b[i]=isapprox(d,0.0,atol=tol_level)
        end
        @test all(b)
    end

    @testset "Direction (s) variations" begin
        face_type =:quad
        close_loop = true

        s=1
        F,V = revolvecurve(Vc,θ;  s=s, n=n,num_steps=num_steps,close_loop=close_loop,face_type=face_type)
        ϕ = [atan(v[2],v[1]) for v in V]
        @test isapprox(minimum(ϕ),0.0,atol=tol_level)
        @test isapprox(maximum(ϕ),π,atol=tol_level)

        s=0
        F,V = revolvecurve(Vc,θ;  s=s, n=n,num_steps=num_steps,close_loop=close_loop,face_type=face_type)
        ϕ = [atan(v[2],v[1]) for v in V]
        @test isapprox(minimum(ϕ),-π/2,atol=tol_level)
        @test isapprox(maximum(ϕ),π/2,atol=tol_level)

        s=-1
        F,V = revolvecurve(Vc,θ;  s=s, n=n,num_steps=num_steps,close_loop=close_loop,face_type=face_type)
        ϕ = [atan(v[2],v[1]) for v in V]
        @test isapprox(minimum(ϕ),-π,atol=tol_level)
        @test isapprox(maximum(ϕ),0.0,atol=tol_level)
    end
    
    @testset "Nothing for num_steps" begin
        face_type =:quad
        close_loop = true
        F,V = revolvecurve(Vc,θ;  s=1, n=n,num_steps=nothing,close_loop=close_loop,face_type=face_type)
        @test length(V)/nc == ceil(Int64,(2*θ)/pointspacingmean(Vc))
    end

    m = [-0.7838430424199713, 0.08108629344330351, -0.6156420208736807]
    @testset "quad" begin        
        face_type =:quad
        close_loop = true
        F,V = revolvecurve(Vc,θ;  s=1, n=m,num_steps=num_steps,close_loop=close_loop,face_type=face_type)

        @test F isa Vector{QuadFace{Int64}}
        @test length(F) == nc*(num_steps-1)
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test isapprox(V[ind], Point3{Float64}[[2.0, 0.0, 0.0], 
        [1.3895382699842485, -0.5610059631967796, 0.06675107120365813], 
        [0.9905120521382602, -1.3374201423606646, 1.1091407653570555], 
        [0.45369147546634303, -0.7152667193399468, 1.2379650904988573], 
        [0.2288198303003941, -0.1271178539023153, 0.965133429366411]], atol=tol_level)

        face_type =:quad
        close_loop = false
        F,V = revolvecurve(Vc,θ;  s=1, n=m,num_steps=num_steps,close_loop=close_loop,face_type=face_type)

        @test F isa Vector{QuadFace{Int64}}
        @test length(F) == (nc-1)*(num_steps-1)
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test isapprox(V[ind], Point3{Float64}[[2.0, 0.0, 0.0], 
        [1.3895382699842485, -0.5610059631967796, 0.06675107120365813], 
        [0.9905120521382602, -1.3374201423606646, 1.1091407653570555], 
        [0.45369147546634303, -0.7152667193399468, 1.2379650904988573], 
        [0.2288198303003941, -0.1271178539023153, 0.965133429366411]], atol=tol_level)
    end

    @testset "tri_slash" begin        
        face_type =:tri_slash
        close_loop = true
        F,V = revolvecurve(Vc,θ;  s=1, n=m,num_steps=num_steps,close_loop=close_loop,face_type=face_type)

        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == nc*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test isapprox(V[ind], Point3{Float64}[[2.0, 0.0, 0.0], 
        [1.3895382699842485, -0.5610059631967796, 0.06675107120365813], 
        [0.9905120521382602, -1.3374201423606646, 1.1091407653570555], 
        [0.45369147546634303, -0.7152667193399468, 1.2379650904988573], 
        [0.2288198303003941, -0.1271178539023153, 0.965133429366411]], atol=tol_level)

        face_type =:tri_slash
        close_loop = false
        F,V = revolvecurve(Vc,θ;  s=1, n=m,num_steps=num_steps,close_loop=close_loop,face_type=face_type)

        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == (nc-1)*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test isapprox(V[ind], Point3{Float64}[[2.0, 0.0, 0.0], 
        [1.3895382699842485, -0.5610059631967796, 0.06675107120365813], 
        [0.9905120521382602, -1.3374201423606646, 1.1091407653570555], 
        [0.45369147546634303, -0.7152667193399468, 1.2379650904988573], 
        [0.2288198303003941, -0.1271178539023153, 0.965133429366411]], atol=tol_level)
    end

    @testset "quad2tri" begin        
        face_type =:quad2tri
        close_loop = true
        F,V = revolvecurve(Vc,θ;  s=1, n=m,num_steps=num_steps,close_loop=close_loop,face_type=face_type)

        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == nc*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test isapprox(V[ind], Point3{Float64}[[2.0, 0.0, 0.0], 
        [1.3895382699842485, -0.5610059631967796, 0.06675107120365813], 
        [0.9905120521382602, -1.3374201423606646, 1.1091407653570555], 
        [0.45369147546634303, -0.7152667193399468, 1.2379650904988573], 
        [0.2288198303003941, -0.1271178539023153, 0.965133429366411]], atol=tol_level)

        face_type =:quad2tri
        close_loop = false
        F,V = revolvecurve(Vc,θ;  s=1, n=m,num_steps=num_steps,close_loop=close_loop,face_type=face_type)

        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == (nc-1)*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test isapprox(V[ind], Point3{Float64}[[2.0, 0.0, 0.0], 
        [1.3895382699842485, -0.5610059631967796, 0.06675107120365813], 
        [0.9905120521382602, -1.3374201423606646, 1.1091407653570555], 
        [0.45369147546634303, -0.7152667193399468, 1.2379650904988573], 
        [0.2288198303003941, -0.1271178539023153, 0.965133429366411]], atol=tol_level)
    end

    @testset "tri" begin        
        face_type =:tri
        close_loop = true
        F,V = revolvecurve(Vc,θ;  s=1, n=m,num_steps=num_steps,close_loop=close_loop,face_type=face_type)

        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == nc*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test isapprox(V[ind], Point3{Float64}[[2.0, 0.0, 0.0], 
        [1.2737434141522277, -0.5142554662637147, 0.06118848193668663], 
        [0.9286050488796189, -1.253831383463123, 1.0398194675222394], 
        [0.45369147546634303, -0.7152667193399468, 1.2379650904988573], 
        [0.2288198303003941, -0.1271178539023153, 0.965133429366411]], atol=tol_level)

        face_type =:tri
        close_loop = false
        F,V = revolvecurve(Vc,θ;  s=1, n=m,num_steps=num_steps,close_loop=close_loop,face_type=face_type)

        @test F isa Vector{TriangleFace{Int64}}
        @test length(F) == (nc-1)*(num_steps-1)*2
        @test V isa Vector{Point3{Float64}}
        @test length(V) == nc*num_steps
        @test isapprox(V[ind], Point3{Float64}[[2.0, 0.0, 0.0], 
        [1.2737434141522277, -0.5142554662637147, 0.06118848193668663], 
        [0.9905120521382602, -1.3374201423606646, 1.1091407653570555], 
        [0.45369147546634303, -0.7152667193399468, 1.2379650904988573], 
        [0.2288198303003941, -0.1271178539023153, 0.965133429366411]], atol=tol_level)
    end

end

@testset "batman" verbose=true begin
    eps_level = 1e-6

    n = 50 # Number of points, first testing an even number
    V = batman(n)
    ind = round.(Int64,range(1,length(V),5))

    @test length(V) == n # Correct length
    @test isa(V,Vector{Point{3,Float64}}) # Correct type
    # Correct coordinates for this case 
    @test isapprox(V[ind],Point3{Float64}[
        [4.959248144765713e-5, -0.6888052278988526, 0.0], 
        [0.9989308524502467, -0.12980802206407838, 0.0], 
        [-0.004310373406159424, 0.3627431045478119, 0.0], 
        [-0.9715858301799283, -0.014868116101481532, 0.0], 
        [-0.022332909569322653, -0.5717894489374953, 0.0]],atol=eps_level)  
          
    V = batman(n; dir=:cw)
    @test length(V) == n # Correct length
    @test isa(V,Vector{Point{3,Float64}})

    V = batman(n; dir=:acw, symmetric=false)
    @test length(V) == n # Correct length
    @test isa(V,Vector{Point{3,Float64}})

    V = batman(n; dir=:acw, symmetric=true)
    @test length(V) == n # Correct length
    @test isa(V,Vector{Point{3,Float64}})

    m = n+1 # force uneven
    V = batman(m; dir=:acw, symmetric=true)
    @test length(V) == m+1 # Correct length
    @test isa(V,Vector{Point{3,Float64}})

    V = batman(m; dir=:cw, symmetric=true)
    @test length(V) == m+1 # Correct length
    @test isa(V,Vector{Point{3,Float64}})
end


@testset "tridisc" verbose=true begin
    eps_level = 1e-6
    r = 2.0
    
    # Test for unrefined hexagon triangulation 
    n = 0 
    F,V = tridisc(r,n)
    @test isa(F,Vector{TriangleFace{Int64}})
    @test isa(V,Vector{Point{3,Float64}})
    @test length(V) == 7
    @test isapprox(maximum(norm.(V)),r,atol=eps_level)

    # Test for defaults    
    F,V = tridisc()
    @test isa(F,Vector{TriangleFace{Int64}})
    @test isa(V,Vector{Point{3,Float64}})
    @test length(V) == 7
    @test isapprox(maximum(norm.(V)),1.0,atol=eps_level)

    # Test for refined hexagon triangulation 
    n = 1 
    F,V = tridisc(r,n)
    @test isa(F,Vector{TriangleFace{Int64}})
    @test isa(V,Vector{Point{3,Float64}})
    @test length(V) == 19
    @test isapprox(maximum(norm.(V)),r,atol=eps_level)
end


@testset "regiontrimesh" verbose=true begin  
    n1 = 120
    r1 = 20.0
    V1 = circlepoints(r1,n1)
    
    n2 = 100
    r2 = 12.0
    V2 = circlepoints(r2,n2)
    V2 = [Point{3,Float64}(v[1]+6,v[2],v[3]) for v in V2]
        
    n3 = 30
    r3 = 4
    V3 = circlepoints(r3,n3)
    V3 = [Point{3,Float64}(v[1]-14,v[2],v[3]) for v in V3]
        
    n4 = 50
    r4 = 7
    V4 = circlepoints(r4,n4)
    V4 = [Point{3,Float64}(v[1]+9,v[2],v[3]) for v in V4]

    n5 = 40
    r5 = 3
    V5 = circlepoints(r5,n5)
    V5 = [Point{3,Float64}(v[1]-2,v[2],v[3]) for v in V5]

    VT = (V1,V2,V3,V4,V5) # Curves
    R = ([1,2,3],[2,4,5],[5]) # Regions 
    P = (1,0.75,0.5)  # Point spacings

    F,V,C = regiontrimesh(VT,R,P)

    @test isa(F,Vector{TriangleFace{Int64}})
    @test isa(V,Vector{Point{3,Float64}})
    @test isa(C,Vector{Float64})
    @test length(F) == length(C) # Color length matches faces

    # Single region case which should trigger treatment of an "on boundary point"
    n = 50
    r = 2.0
    V = circlepoints(r,n;dir=:acw) 
    pointSpacing = pointspacingmean(V)

    VT = (V,)
    R = ([1],)
    P = (pointSpacing)
    F,V,C = regiontrimesh(VT,R,P)
    @test isa(F,Vector{TriangleFace{Int64}})
    @test isa(V,Vector{Point{3,Float64}})
    @test isa(C,Vector{Float64})
    @test length(F) == length(C)

end


@testset "scalesimplex" verbose=true begin
    eps_level = 1e-6

    M = cube(sqrt(3))
    Fq = faces(M)
    Vq = coordinates(M)

    M = tetrahedron(1.0)
    Ft = faces(M)
    Vt = coordinates(M)

    @testset "Single scaling value applied to quads" begin
        Fs,Vs = scalesimplex(Fq,Vq,0.5)
        @test typeof(Fq) == typeof(Fs)
        @test typeof(Vq) == typeof(Vs)
        @test length(Fq) == length(Fs)
        @test length(Vs) == length(Fq)*length(Fq[1]) 
        @test isapprox(Vs,Point{3, Float64}[[-0.5, -0.5, -1.0], [-0.5, 0.5, -1.0], 
        [0.5, 0.5, -1.0], [0.5, -0.5, -1.0], [0.5, -0.5, 1.0], [0.5, 0.5, 1.0], 
        [-0.5, 0.5, 1.0], [-0.5, -0.5, 1.0], [-1.0, -0.5, 0.5], [-1.0, 0.5, 0.5], 
        [-1.0, 0.5, -0.5], [-1.0, -0.5, -0.5], [-0.5, 1.0, 0.5], [0.5, 1.0, 0.5], 
        [0.5, 1.0, -0.5], [-0.5, 1.0, -0.5], [1.0, 0.5, 0.5], [1.0, -0.5, 0.5], 
        [1.0, -0.5, -0.5], [1.0, 0.5, -0.5], [0.5, -1.0, 0.5], [-0.5, -1.0, 0.5], 
        [-0.5, -1.0, -0.5], [0.5, -1.0, -0.5]],atol=eps_level)
    end

    @testset "Per face scaling value applied to quads" begin
        Fs,Vs = scalesimplex(Fq,Vq,0.5*ones(length(Fq)))
        @test typeof(Fq) == typeof(Fs)
        @test typeof(Vq) == typeof(Vs)
        @test length(Fq) == length(Fs)
        @test length(Vs) == length(Fq)*length(Fq[1]) 
        @test isapprox(Vs,Point{3, Float64}[[-0.5, -0.5, -1.0], [-0.5, 0.5, -1.0], 
        [0.5, 0.5, -1.0], [0.5, -0.5, -1.0], [0.5, -0.5, 1.0], [0.5, 0.5, 1.0], 
        [-0.5, 0.5, 1.0], [-0.5, -0.5, 1.0], [-1.0, -0.5, 0.5], [-1.0, 0.5, 0.5], 
        [-1.0, 0.5, -0.5], [-1.0, -0.5, -0.5], [-0.5, 1.0, 0.5], [0.5, 1.0, 0.5], 
        [0.5, 1.0, -0.5], [-0.5, 1.0, -0.5], [1.0, 0.5, 0.5], [1.0, -0.5, 0.5], 
        [1.0, -0.5, -0.5], [1.0, 0.5, -0.5], [0.5, -1.0, 0.5], [-0.5, -1.0, 0.5], 
        [-0.5, -1.0, -0.5], [0.5, -1.0, -0.5]],atol=eps_level)
    end

    @testset "Per face scaling value applied to quads" begin
        Fs,Vs = scalesimplex(Fq,Vq,0.5*ones(length(Vq)))
        @test typeof(Fq) == typeof(Fs)
        @test typeof(Vq) == typeof(Vs)
        @test length(Fq) == length(Fs)
        @test length(Vs) == length(Fq)*length(Fq[1]) 
        @test isapprox(Vs,Point{3, Float64}[[-0.5, -0.5, -1.0], [-0.5, 0.5, -1.0], 
        [0.5, 0.5, -1.0], [0.5, -0.5, -1.0], [0.5, -0.5, 1.0], [0.5, 0.5, 1.0], 
        [-0.5, 0.5, 1.0], [-0.5, -0.5, 1.0], [-1.0, -0.5, 0.5], [-1.0, 0.5, 0.5], 
        [-1.0, 0.5, -0.5], [-1.0, -0.5, -0.5], [-0.5, 1.0, 0.5], [0.5, 1.0, 0.5], 
        [0.5, 1.0, -0.5], [-0.5, 1.0, -0.5], [1.0, 0.5, 0.5], [1.0, -0.5, 0.5], 
        [1.0, -0.5, -0.5], [1.0, 0.5, -0.5], [0.5, -1.0, 0.5], [-0.5, -1.0, 0.5], 
        [-0.5, -1.0, -0.5], [0.5, -1.0, -0.5]],atol=eps_level)
    end

    @testset "Single scaling value applied to triangles" begin
        Fs,Vs = scalesimplex(Ft,Vt,0.5)
        @test typeof(Ft) == typeof(Fs)
        @test typeof(Vt) == typeof(Vs)
        @test length(Ft) == length(Fs)
        @test length(Vs) == length(Ft)*length(Ft[1]) 
        @test isapprox(Vs, Point{3, Float64}[[-0.4082482904638631, -0.3928371006591931, -0.11111111111111109], 
        [0.4082482904638631, -0.3928371006591931, -0.11111111111111109], 
        [0.0, -0.15713484026367724, 0.5555555555555556], 
        [0.0, 0.47140452079103173, -0.3333333333333333], 
        [0.4082482904638631, -0.23570226039551587, -0.3333333333333333], 
        [-0.4082482904638631, -0.23570226039551587, -0.3333333333333333], 
        [0.13608276348795437, 0.5499719409228703, -0.11111111111111109], 
        [0.13608276348795437, 0.07856742013183862, 0.5555555555555556], 
        [0.5443310539518174, -0.15713484026367724, -0.11111111111111109], 
        [-0.13608276348795437, 0.5499719409228703, -0.11111111111111109], 
        [-0.5443310539518174, -0.15713484026367724, -0.11111111111111109], 
        [-0.13608276348795437, 0.07856742013183862, 0.5555555555555556]],atol=eps_level)
    end

end


@testset "subcurve" verbose=true begin
    eps_level = 1e-6

    t = range(0.0,2*π,5) 
    V = [Point{3,Float32}(ti,cos(ti),0.0) for ti in t] # Data values

    n = 1
    Vn = subcurve(V,n)
    @test typeof(V) == typeof(Vn)
    @test length(Vn) == length(V) + (length(V)-1)*n 

    t = range(0.0,2*π,5) 
    V = [Point{3,Float64}(ti,cos(ti),0.0) for ti in t] # Data values

    n = 2
    Vn = subcurve(V,n)
    @test typeof(V) == typeof(Vn)
    @test length(Vn) == length(V) + (length(V)-1)*n 

    @test isapprox(Vn,Point{3, Float64}[[0.0, 1.0, 0.0], [0.5235987755982988, 0.6666666666666667, 0.0], 
    [1.0471975511965976, 0.3333333333333334, 0.0], [1.5707963267948966, 6.123233995736766e-17, 0.0], 
    [2.0943951023931953, -0.33333333333333326, 0.0], [2.617993877991494, -0.6666666666666666, 0.0], 
    [3.141592653589793, -1.0, 0.0], [3.6651914291880923, -0.6666666666666667, 0.0], 
    [4.188790204786391, -0.3333333333333335, 0.0], [4.71238898038469, -1.8369701987210297e-16, 0.0], 
    [5.235987755982988, 0.33333333333333315, 0.0], [5.759586531581287, 0.6666666666666665, 0.0], 
    [6.283185307179586, 1.0, 0.0]],atol=eps_level)
end

@testset "dualclad" verbose=true begin
    eps_level = 1e-6

    @testset "Triangles, single scale" begin
        M = tetrahedron(1.0)
        F = faces(M)
        V = coordinates(M)
        E = meshedges(F; unique_only = true)
        Fs,Fsq,Vs = dualclad(F,V,0.5; connectivity=:face)
        @test length(Fs) == length(F)
        @test length(Fsq) == length(E)
        @test length(Vs) == length(F)*length(F[1])
        @test isapprox(Vs,Point{3, Float64}[[-0.4082482904638631, -0.3928371006591931, -0.11111111111111109], 
        [0.4082482904638631, -0.3928371006591931, -0.11111111111111109], 
        [0.0, -0.15713484026367724, 0.5555555555555556], [0.0, 0.47140452079103173, -0.3333333333333333], 
        [0.4082482904638631, -0.23570226039551587, -0.3333333333333333], 
        [-0.4082482904638631, -0.23570226039551587, -0.3333333333333333], 
        [0.13608276348795437, 0.5499719409228703, -0.11111111111111109], 
        [0.13608276348795437, 0.07856742013183862, 0.5555555555555556], 
        [0.5443310539518174, -0.15713484026367724, -0.11111111111111109],
         [-0.13608276348795437, 0.5499719409228703, -0.11111111111111109], 
         [-0.5443310539518174, -0.15713484026367724, -0.11111111111111109], 
         [-0.13608276348795437, 0.07856742013183862, 0.5555555555555556]],atol=eps_level)

        Fs,Fsq,Vs = dualclad(F,V,0.5; connectivity=:edge)
        @test length(Fs) == length(F)
        @test length(Fsq) == length(F)*length(F[1])
        @test length(Vs) == length(F)*length(F[1]) + length(E)*2
        @test isapprox(Vs,Point{3, Float64}[[-0.4082482904638631, -0.3928371006591931, -0.11111111111111109], 
        [0.4082482904638631, -0.3928371006591931, -0.11111111111111109], 
        [0.0, -0.15713484026367724, 0.5555555555555556], 
        [0.0, 0.47140452079103173, -0.3333333333333333], 
        [0.4082482904638631, -0.23570226039551587, -0.3333333333333333], 
        [-0.4082482904638631, -0.23570226039551587, -0.3333333333333333], 
        [0.13608276348795437, 0.5499719409228703, -0.11111111111111109], 
        [0.13608276348795437, 0.07856742013183862, 0.5555555555555556], 
        [0.5443310539518174, -0.15713484026367724, -0.11111111111111109], 
        [-0.13608276348795437, 0.5499719409228703, -0.11111111111111109], 
        [-0.5443310539518174, -0.15713484026367724, -0.11111111111111109], 
        [-0.13608276348795437, 0.07856742013183862, 0.5555555555555556], 
        [-0.4082482904638631, -0.47140452079103173, -0.3333333333333333], 
        [0.4082482904638631, -0.47140452079103173, -0.3333333333333333], 
        [0.20412414523193154, 0.5892556509887896, -0.3333333333333333], 
        [0.6123724356957946, -0.11785113019775792, -0.3333333333333333], 
        [0.0, 0.7071067811865476, 0.0], [0.0, 0.23570226039551587, 0.6666666666666667], 
        [-0.20412414523193154, 0.5892556509887896, -0.3333333333333333], 
        [-0.6123724356957946, -0.11785113019775792, -0.3333333333333333], 
        [0.6123724356957946, -0.3535533905932738, 0.0], 
        [0.20412414523193154, -0.11785113019775793, 0.6666666666666667], 
        [-0.6123724356957946, -0.3535533905932738, 0.0], 
        [-0.20412414523193154, -0.11785113019775793, 0.6666666666666667]],atol=eps_level)
    end    

    @testset "Triangles, per face and per edge scale" begin
        M = tetrahedron(1.0)
        F = faces(M)
        V = coordinates(M)
        E = meshedges(F; unique_only = true)
        Fs,Fsq,Vs = dualclad(F,V,0.5*ones(length(F)); connectivity=:face)
        @test length(Fs) == length(F)
        @test length(Fsq) == length(E)
        @test length(Vs) == length(F)*length(F[1])
        @test isapprox(Vs,Point{3, Float64}[[-0.4082482904638631, -0.3928371006591931, -0.11111111111111109], 
        [0.4082482904638631, -0.3928371006591931, -0.11111111111111109], 
        [0.0, -0.15713484026367724, 0.5555555555555556], [0.0, 0.47140452079103173, -0.3333333333333333], 
        [0.4082482904638631, -0.23570226039551587, -0.3333333333333333], 
        [-0.4082482904638631, -0.23570226039551587, -0.3333333333333333], 
        [0.13608276348795437, 0.5499719409228703, -0.11111111111111109], 
        [0.13608276348795437, 0.07856742013183862, 0.5555555555555556], 
        [0.5443310539518174, -0.15713484026367724, -0.11111111111111109],
         [-0.13608276348795437, 0.5499719409228703, -0.11111111111111109], 
         [-0.5443310539518174, -0.15713484026367724, -0.11111111111111109], 
         [-0.13608276348795437, 0.07856742013183862, 0.5555555555555556]],atol=eps_level)

        Fs,Fsq,Vs = dualclad(F,V,0.5*ones(length(V)); connectivity=:edge)
        @test length(Fs) == length(F)
        @test length(Fsq) == length(F)*length(F[1])
        @test length(Vs) == length(F)*length(F[1]) + length(E)*2
        @test isapprox(Vs,Point{3, Float64}[[-0.4082482904638631, -0.3928371006591931, -0.11111111111111109], 
        [0.4082482904638631, -0.3928371006591931, -0.11111111111111109], 
        [0.0, -0.15713484026367724, 0.5555555555555556], 
        [0.0, 0.47140452079103173, -0.3333333333333333], 
        [0.4082482904638631, -0.23570226039551587, -0.3333333333333333], 
        [-0.4082482904638631, -0.23570226039551587, -0.3333333333333333], 
        [0.13608276348795437, 0.5499719409228703, -0.11111111111111109], 
        [0.13608276348795437, 0.07856742013183862, 0.5555555555555556], 
        [0.5443310539518174, -0.15713484026367724, -0.11111111111111109], 
        [-0.13608276348795437, 0.5499719409228703, -0.11111111111111109], 
        [-0.5443310539518174, -0.15713484026367724, -0.11111111111111109], 
        [-0.13608276348795437, 0.07856742013183862, 0.5555555555555556], 
        [-0.4082482904638631, -0.47140452079103173, -0.3333333333333333], 
        [0.4082482904638631, -0.47140452079103173, -0.3333333333333333], 
        [0.20412414523193154, 0.5892556509887896, -0.3333333333333333], 
        [0.6123724356957946, -0.11785113019775792, -0.3333333333333333], 
        [0.0, 0.7071067811865476, 0.0], [0.0, 0.23570226039551587, 0.6666666666666667], 
        [-0.20412414523193154, 0.5892556509887896, -0.3333333333333333], 
        [-0.6123724356957946, -0.11785113019775792, -0.3333333333333333], 
        [0.6123724356957946, -0.3535533905932738, 0.0], 
        [0.20412414523193154, -0.11785113019775793, 0.6666666666666667], 
        [-0.6123724356957946, -0.3535533905932738, 0.0], 
        [-0.20412414523193154, -0.11785113019775793, 0.6666666666666667]],atol=eps_level)
    end    
    
    @testset "Quads, single scale" begin
        M = cube(√3)
        F = faces(M)
        V = coordinates(M)
        E = meshedges(F; unique_only = true)
        Fs,Fsq,Vs = dualclad(F,V,0.5; connectivity=:face)
        @test length(Fs) == length(F)
        @test length(Fsq) == length(E)
        @test length(Vs) == length(F)*length(F[1])
        @test isapprox(Vs,Point{3, Float64}[[-0.5, -0.5, -1.0], [-0.5, 0.5, -1.0],
         [0.5, 0.5, -1.0], [0.5, -0.5, -1.0], [0.5, -0.5, 1.0], [0.5, 0.5, 1.0], 
         [-0.5, 0.5, 1.0], [-0.5, -0.5, 1.0], [-1.0, -0.5, 0.5], [-1.0, 0.5, 0.5], 
         [-1.0, 0.5, -0.5], [-1.0, -0.5, -0.5], [-0.5, 1.0, 0.5], [0.5, 1.0, 0.5],
          [0.5, 1.0, -0.5], [-0.5, 1.0, -0.5], [1.0, 0.5, 0.5], [1.0, -0.5, 0.5],
           [1.0, -0.5, -0.5], [1.0, 0.5, -0.5], [0.5, -1.0, 0.5], [-0.5, -1.0, 0.5], 
           [-0.5, -1.0, -0.5], [0.5, -1.0, -0.5]],atol=eps_level)

        Fs,Fsq,Vs = dualclad(F,V,0.5; connectivity=:edge)
        @test length(Fs) == length(F)
        @test length(Fsq) == length(F)*length(F[1])
        @test length(Vs) == length(F)*length(F[1]) + length(E)*2
        @test isapprox(Vs,Point{3, Float64}[[-0.5, -0.5, -1.0], [-0.5, 0.5, -1.0], 
        [0.5, 0.5, -1.0], [0.5, -0.5, -1.0], [0.5, -0.5, 1.0], [0.5, 0.5, 1.0], 
        [-0.5, 0.5, 1.0], [-0.5, -0.5, 1.0], [-1.0, -0.5, 0.5], [-1.0, 0.5, 0.5], 
        [-1.0, 0.5, -0.5], [-1.0, -0.5, -0.5], [-0.5, 1.0, 0.5], [0.5, 1.0, 0.5], 
        [0.5, 1.0, -0.5], [-0.5, 1.0, -0.5], [1.0, 0.5, 0.5], [1.0, -0.5, 0.5], 
        [1.0, -0.5, -0.5], [1.0, 0.5, -0.5], [0.5, -1.0, 0.5], [-0.5, -1.0, 0.5], 
        [-0.5, -1.0, -0.5], [0.5, -1.0, -0.5], [-1.0, -0.5, -1.0], [-1.0, 0.5, -1.0], 
        [1.0, -0.5, 1.0], [1.0, 0.5, 1.0], [-1.0, -0.5, 1.0], [-1.0, 0.5, 1.0], 
        [-0.5, 1.0, 1.0], [0.5, 1.0, 1.0], [0.5, -1.0, 1.0], [-0.5, -1.0, 1.0], 
        [-0.5, 1.0, -1.0], [0.5, 1.0, -1.0], [-1.0, 1.0, 0.5], [-1.0, 1.0, -0.5], 
        [1.0, 1.0, 0.5], [1.0, 1.0, -0.5], [1.0, -1.0, 0.5], [1.0, -1.0, -0.5], 
        [-1.0, -1.0, 0.5], [-1.0, -1.0, -0.5], [1.0, 0.5, -1.0], [1.0, -0.5, -1.0],
         [-0.5, -1.0, -1.0], [0.5, -1.0, -1.0]],atol=eps_level)
    end    
end