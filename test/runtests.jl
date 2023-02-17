using SyntheticEddyMethod
using Test
using Revise
@testset "SyntheticEddyMethod.jl" begin
    @test SyntheticEddyMethod.greet_your_package_name() == "Hello YourPackageName!"
    @test SyntheticEddyMethod.greet_your_package_name() != "Hello world!"
end
