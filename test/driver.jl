
using Statistics
function driver_test(TI::Float64; Nt = 1000, s_fun = tent_fun, DFSEM = false)
σ = 0.1 #eddy dimensions, the same in all the directions
b = 5.0
a = 0.0

#Defining the Virtual Box domain
x = -σ:0.1:+σ 
y = collect(a:0.1:b)
z = collect(a:0.1:b)


Vboxinfo = VirtualBox(y,z,σ; shape_fun = s_fun)

@test Vboxinfo.σ == [σ,σ,σ]
@test Vboxinfo.X_end == σ
@test Vboxinfo.X_start == -σ
@test Vboxinfo.Y_end == b+σ
@test Vboxinfo.Y_start == a-σ
@test Vboxinfo.Z_end == b+σ
@test Vboxinfo.Z_start == a-σ

Vboxinfo.N = 100
@test Vboxinfo.N == 100

Vboxinfo = VirtualBox(y,z,σ; shape_fun = s_fun)
N = Vboxinfo.N

dt = 0.01

U₀ = 1.0 #Convective Velocity
#TI #turbulence intensity


#Isotropic turbulence
u_p = (U₀ * TI)^2

Re_stress = [u_p 0.0 0.0; 
            0.0 u_p 0.0;
            0.0 0.0 u_p]


Eddies = initialize_eddies(Vboxinfo)

# Or you can use the wrapper function
Re_stress_2, Eddies = initialize_eddies(U₀, TI, Vboxinfo)

@test Re_stress == Re_stress_2

#Computing the velocity in the middle of the VirtualBox domain
vector_points = create_vector_points(x, y, z)

@test length(vector_points) == length(x)*length(y)*length(z)


#Computing the velocity in the middle of the VirtualBox domain
vector_points = [[0.0, b/2, b/2]]

#Defining how many time interval

Nt = Nt
U = zeros(Nt, 3)

u_f = compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo, Re_stress; DFSEM = DFSEM)
@test typeof(u_f) == Vector{Vector{Float64}}
@test length(u_f) == length(vector_points)

for i = 1:1:Nt
    u_f = compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo, Re_stress; DFSEM = DFSEM)
    U[i,:] = u_f[1] #Its a vector of vector
end

#The deviation standard should approach the turbulence intensity
    s1 = Statistics.std(U[:,1])
    s2 = Statistics.std(U[:,2])
    s3 = Statistics.std(U[:,3])

return mean([s1,s2,s3])

end


function utilities_test()
    return isapprox(compute_Ek([[1.2, 0.1, 0.3]], 1.0)[1], (0.2^2 + 0.1^2 + 0.3^2 )*0.5)
end