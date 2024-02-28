using Test, Comodo 

@testset "dist" verbose = true begin 

    v1 = Float64[0, 0, 0]
    v2 = Float64[0, 0, 5]

    result = dist(v1, v2)

    @test result == [0.0 0.0 5.0; 0.0 0.0 5.0; 0.0 0.0 5.0]
end 