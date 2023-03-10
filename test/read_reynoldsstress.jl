using XLSX
 function read_re_test()
   ftest = ["Re_ch.xlsx", "Re_stress.xlsx"]
   dimstest = [(:Z, :Y), (:Y,)]
   for (fname,dims) in zip(ftest,dimstest)
    reynolds_stress_file = joinpath(@__DIR__,"Data",fname)
    Re_from_file = get_reynolds_stress_from_file(reynolds_stress_file; dims = dims)
    @test typeof(Re_from_file) == Reynolds_stress_interpolator

    σ = 0.1 #eddy dimensions, the same in all the directions
    a = 0.2

    #Defining the Virtual Box domain
    y = collect(a:0.1:1.8)
    z = collect(a:0.1:6)


    Vboxinfo = VirtualBox(y,z,σ)
    dt = 0.01
    U₀ = 1.0 #Convective Velocity
    Eddies = initialize_eddies(Vboxinfo)
    vector_points = [[0.0, 0.2, 3.0], [0.0, 0.9, 4.0]]
    
    u_f = compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo, Re_from_file)

    @test typeof(u_f) == Vector{Vector{Float64}}

   end
 end

