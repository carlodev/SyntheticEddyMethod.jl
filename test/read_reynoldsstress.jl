using XLSX
function read_re_test()
    reynolds_stress_file = joinpath(@__DIR__,"Data","Re_ch.xlsx")
    A_from_file = get_reynolds_stress_from_file(reynolds_stress_file)
    @test typeof(A_from_file) == Reynolds_stress_interpolator
end
