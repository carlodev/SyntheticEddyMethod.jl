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


###
### Verify that the Reynolds stress created is close to the provide one
# fname = "Re_stress.xlsx"
# dims = (:Y,)
# reynolds_stress_file = joinpath(@__DIR__,"Data",fname)
# Re_from_file = get_reynolds_stress_from_file(reynolds_stress_file; dims = dims)

# σ = 0.05 #eddy dimensions, the same in all the directions

# using SyntheticEddyMethod

# #Defining the Virtual Box domain
# y = collect(-1:0.1:1)
# z = collect(-2.1:0.1:2.1)

# Vboxinfo = VirtualBox(y,z,σ)
# dt = 0.01
# U₀ = 1.0 #Convective Velocity
# Eddies = initialize_eddies(Vboxinfo)
# point = [0.0, 0.2, 0.0]

# Re = [2.45e-1  -7.03e-2 2.68e-3;
#       -7.03e-2 2.3e-1 -4.46e-5;
#       2.68e-3  -4.46e-5 4.39e-1]

# Re = [2.0  0.001 -0.001;
#       0.001 1.0 0.001;
#       -0.001  0.001 2.0]

# eig_val, eig_vec = eigen(Re)
# sum(eig_val) - 2*maximum(eig_val)

# Nt = 10000
# U = zeros(Nt, 3)
# time_vec = collect(0:dt:dt*(Nt-1))
# for i = 1:1:Nt
#   u_f = compute_fluct(point, time_vec[i], Eddies, U₀, Vboxinfo, Re; DFSEM = false)
#   U[i,:] = u_f #A vector of vector
# end

# std(U[:,1])^2
# std(U[:,2])^2
# std(U[:,3])^2

# Statistics.mean((U[:,1] .-U₀) .* U[:,2])
# Statistics.mean((U[:,1] .-U₀) .* U[:,3])
# Statistics.mean(U[:,2] .* U[:,3])
