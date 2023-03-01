using SyntheticEddyMethod
using Test
using Plots
@testset "SyntheticEddyMethod.jl" begin
    include("driver.jl")
    @test isapprox(driver_test(0.1), 0.1; rtol =0.2)
    @test isapprox(driver_test(0.01), 0.01; rtol =0.2)
    @test isapprox(driver_test(0.05), 0.05; rtol =0.2)

    @test isapprox(driver_test(0.01; s_fun = step_fun), 0.01; rtol =0.2)
    @test isapprox(driver_test(0.01; s_fun = trunc_gauss_fun), 0.01; rtol =0.2)


end    
