using SyntheticEddyMethod
using Test
@testset "SyntheticEddyMethod.jl" begin
    include("driver.jl")
    @test isapprox(driver_test(0.1), 0.1; rtol =0.2)
    @test isapprox(driver_test(0.01), 0.01; rtol =0.2)
    @test isapprox(driver_test(0.05), 0.05; rtol =0.2)
end    
