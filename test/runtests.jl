using Test, Comodo, GeometryBasics

@testset "dist" verbose = true begin

	v1 = Float64[0, 0, 0]
	v2 = Float64[0, 0, 5]

	result = dist(v1, v2)

	@test result == [0.0 0.0 5.0; 0.0 0.0 5.0; 0.0 0.0 5.0]
end

@testset "gridpoints" verbose = true begin

	@testset "with 1 vector" begin
		a = Float64[1, 2, 3]

		expected = Point3{Float64}[
			[1.0, 1.0, 1.0],
			[1.0, 1.0, 2.0],
			[1.0, 1.0, 3.0],
			[1.0, 2.0, 1.0],
			[1.0, 2.0, 2.0],
			[1.0, 2.0, 3.0],
			[1.0, 3.0, 1.0],
			[1.0, 3.0, 2.0],
			[1.0, 3.0, 3.0],
			[2.0, 1.0, 1.0],
			[2.0, 1.0, 2.0],
			[2.0, 1.0, 3.0],
			[2.0, 2.0, 1.0],
			[2.0, 2.0, 2.0],
			[2.0, 2.0, 3.0],
			[2.0, 3.0, 1.0],
			[2.0, 3.0, 2.0],
			[2.0, 3.0, 3.0],
			[3.0, 1.0, 1.0],
			[3.0, 1.0, 2.0],
			[3.0, 1.0, 3.0],
			[3.0, 2.0, 1.0],
			[3.0, 2.0, 2.0],
			[3.0, 2.0, 3.0],
			[3.0, 3.0, 1.0],
			[3.0, 3.0, 2.0],
			[3.0, 3.0, 3.0],
		]

		result = gridpoints(a)

		@test result == expected
	end

	@testset "with 2 vectors" begin
		a = Float64[1, 2, 3]
		b = Float64[2, 3, 5]

		expected = Point3{Float64}[
			[1.0, 2.0, 1.0],
			[1.0, 2.0, 2.0],
			[1.0, 2.0, 3.0],
			[1.0, 3.0, 1.0],
			[1.0, 3.0, 2.0],
			[1.0, 3.0, 3.0],
			[1.0, 5.0, 1.0],
			[1.0, 5.0, 2.0],
			[1.0, 5.0, 3.0],
			[2.0, 2.0, 1.0],
			[2.0, 2.0, 2.0],
			[2.0, 2.0, 3.0],
			[2.0, 3.0, 1.0],
			[2.0, 3.0, 2.0],
			[2.0, 3.0, 3.0],
			[2.0, 5.0, 1.0],
			[2.0, 5.0, 2.0],
			[2.0, 5.0, 3.0],
			[3.0, 2.0, 1.0],
			[3.0, 2.0, 2.0],
			[3.0, 2.0, 3.0],
			[3.0, 3.0, 1.0],
			[3.0, 3.0, 2.0],
			[3.0, 3.0, 3.0],
			[3.0, 5.0, 1.0],
			[3.0, 5.0, 2.0],
			[3.0, 5.0, 3.0],
		]

		result = gridpoints(a, b)

		@test result == expected
	end

	@testset "with 3 vectors" begin
		a = Float64[1, 2, 3]
		b = Float64[2, 3, 5]
		c = Float64[5, 6, 4]

		expected = GeometryBasics.Point3{Float64}[
			[1.0, 2.0, 5.0],
			[1.0, 2.0, 6.0],
			[1.0, 2.0, 4.0],
			[1.0, 3.0, 5.0],
			[1.0, 3.0, 6.0],
			[1.0, 3.0, 4.0],
			[1.0, 5.0, 5.0],
			[1.0, 5.0, 6.0],
			[1.0, 5.0, 4.0],
			[2.0, 2.0, 5.0],
			[2.0, 2.0, 6.0],
			[2.0, 2.0, 4.0],
			[2.0, 3.0, 5.0],
			[2.0, 3.0, 6.0],
			[2.0, 3.0, 4.0],
			[2.0, 5.0, 5.0],
			[2.0, 5.0, 6.0],
			[2.0, 5.0, 4.0],
			[3.0, 2.0, 5.0],
			[3.0, 2.0, 6.0],
			[3.0, 2.0, 4.0],
			[3.0, 3.0, 5.0],
			[3.0, 3.0, 6.0],
			[3.0, 3.0, 4.0],
			[3.0, 5.0, 5.0],
			[3.0, 5.0, 6.0],
			[3.0, 5.0, 4.0],
		]

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


@testset "geosphere(3,1.0)" begin
	F, V = geosphere(3, 1.0)

	@test F isa Vector{TriangleFace{Int64}}
	@test length(F) == 1280

	@test V isa Vector{Point3{Float64}}
	@test length(V) == 642
end

@testset "simplexcenter" begin

	eps = 0.001

	F, V = geosphere(3, 1.0)
	VC = simplexcenter(F, V)

	@test typeof(VC) == Vector{Point3{Float64}}
	@test length(VC) == 1280
	@test isapprox(VC[1], [-0.35520626817942325, 0.0, -0.9299420831107401], atol = eps)

end


@testset "remove unused vertices" begin
	F, V = geosphere(3, 1.0)
	VC = simplexcenter(F, V)
	F = [F[i] for i in findall(map(v -> v[3] > 0, VC))] # Remove some faces
	F, V = remove_unused_vertices(F, V)

	@test F isa Vector{TriangleFace{Int64}}
	@test length(F) == 624

	@test V isa Vector{Point3{Float64}}
	@test length(V) == 337
end


@testset "boundry edges" begin
	F, V = geosphere(3, 1.0)
	VC = simplexcenter(F, V)
	F = [F[i] for i in findall(map(v -> v[3] > 0, VC))] # Remove some faces
	F, V = remove_unused_vertices(F, V)
	Eb = boundaryedges(F)

	@test typeof(Eb) == Vector{LineFace{Int64}}
	@test length(Eb) == 48
	@test Eb[1] == [272, 205]
end

@testset "edges to curve" begin
	F, V = geosphere(3, 1.0)
	VC = simplexcenter(F, V)
	F = [F[i] for i in findall(map(v -> v[3] > 0, VC))] # Remove some faces
	F, V = remove_unused_vertices(F, V)
	Eb = boundaryedges(F)
	ind = edges2curve(Eb)

	@test length(ind) == 49
	@test ind ==
		  [272, 205, 329, 168, 309, 184, 270, 145, 306, 220, 334, 223, 305, 138, 320, 204, 292, 232, 336, 234, 311, 203, 269, 115, 303, 194, 321, 133, 312, 207, 271, 164, 308, 209, 330, 231, 307, 154, 327, 206, 301, 240, 335, 229, 310, 196, 304, 208, 272]
end


@testset "mesh" begin
	F, V = geosphere(3, 1.0)
	VC = simplexcenter(F, V)
	F = [F[i] for i in findall(map(v -> v[3] > 0, VC))] # Remove some faces
	F, V = remove_unused_vertices(F, V)

	M = GeometryBasics.Mesh(V, F)

	@test typeof(M) == Mesh{3, Float64, GeometryBasics.Ngon{3, Float64, 3, Point3{Float64}}, SimpleFaceView{3, Float64, 3, Int64, Point3{Float64}, TriangleFace{Int64}}}
	@test length(M) == 624
	@test M[1][1] == [-0.5642542117657715, -0.5133754412304479, 0.6465777917977317]

end


@testset "circle points" begin 

    @testset "with value" begin
        V1 = circlepoints(1.0, 40)

        @test V1 isa Vector{Point3{Float64}}
        @test length(V1) == 40
        @test V1[1] == [1.0, 0.0, 0.0]
    end 

    @testset "with function" begin 
        r = 1.0
        n = 40
        rFun(t) = r + 0.5.*sin(3*t)
        V2 = circlepoints(rFun, n)
        @test V2 isa Vector{Point3{Float64}}
        @test length(V2) == 40
        @test V2[1] == [1.0, 0.0, 0.0]
    end 

end 