using SyntheticEddyMethod
using Test
using Plots



@testset "VirtualBox" begin include("TestDrivers/driver_virtualbox.jl") end

@testset "Fluctuations" begin include("TestDrivers/driver_fluctuations.jl") 
    fluctuations_test(0.1)
    fluctuations_test(0.01)
    fluctuations_test(0.05)
    fluctuations_test(0.01; s_fun = step_fun)
    fluctuations_test(0.01; s_fun = trunc_gauss_fun)
    fluctuations_test(0.01; s_fun = DFSEM_fun, DFSEM = true)
    test_div_null(0.001)
end

@testset "ReynoldsStress" begin include("TestDrivers/driver_reynoldsstress.jl") 
    read_re_test()
    test_anisotropic_reynolds()
end 


@testset "Utilities" begin include("TestDrivers/driver_utilities.jl") end