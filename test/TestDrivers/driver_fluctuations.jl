
using Statistics
using LinearAlgebra



function fluctuations_test(TI::Float64; Nt = 1000, s_fun = tent_fun, DFSEM = false)
σ = 0.1 #eddy dimensions, the same in all the directions
b = 5.0
a = 0.0

#Defining the Virtual Box domain
x = -σ:0.1:+σ 
y = collect(a:0.1:b)
z = collect(a:0.1:b)


Vboxinfo = VirtualBox(y,z,σ; shape_fun = s_fun)

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
point = [0.0, b/2, b/2]

#Defining how many time interval

U = zeros(Nt, 3)

u_f = compute_fluct(point, dt, Eddies, U₀, Vboxinfo, Re_stress; DFSEM = DFSEM)
@test typeof(u_f) == Vector{Float64}
@test length(u_f) == 3

time_vec = collect(0:dt:dt*(Nt-1))
@info "Convecting Eddies for $Nt time steps"
for i = 1:1:Nt
    u_f = compute_fluct(point, time_vec[i], Eddies, U₀, Vboxinfo, Re_stress; DFSEM = DFSEM)
    U[i,:] = u_f #A vector of vector
end

#The deviation standard should approach the turbulence intensity
    s1 = Statistics.std(U[:,1])
    s2 = Statistics.std(U[:,2])
    s3 = Statistics.std(U[:,3])

TI_computed = mean([s1,s2,s3])

#TEST TI computed == TI provided by the user
@test isapprox(TI_computed, TI; rtol =0.2)
 

end

function test_div_null(tol=0.1)
    σ = 0.1
    b = 5.0
    a = 0.0
    y = collect(a:0.1:b)
    z = collect(a:0.1:b)
    
    Vboxinfo = VirtualBox(y,z,σ)
    dt = 0.001
    U₀ = 1.0
    TI = 0.01 #turbulence intensity
    norm_div = 1
    iter = 0
    iter_max = 100
    while (norm_div > tol || isnan(norm_div)) && iter<iter_max
        Re_stress, Eddies = initialize_eddies(U₀, TI, Vboxinfo)
        
        dl = 0.0001
        vector_points = [[0.0, b / 2, b / 2], [dl, b / 2, b / 2], [0.0, b / 2 + dl, b / 2], [0.0, b / 2, b / 2 + dl]]
        u_fluct = compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo, Re_stress; DFSEM = true)

        dudx = (u_fluct[2][1] - u_fluct[1][1]) / dl
        dvdy = (u_fluct[3][2] - u_fluct[1][2]) / dl
        dwdz = (u_fluct[4][3] - u_fluct[1][3]) / dl
        grad_norm = norm([dudx, dvdy, dwdz])
        div_val = dudx + dvdy + dwdz
        norm_div = abs(div_val/grad_norm)
        iter = iter+1
    
    end
    @info "iter = $iter for satifying the requirement $tol"
    @test norm_div<tol

end


