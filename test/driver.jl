
using Statistics
function driver_test(TI::Float64; Nt = 2000)
σ = 0.1 #eddy dimensions, the same in all the directions
b = 5.0
a = 0.0

#Defining the Virtual Box domain
x = -σ:0.1:+σ 
y = collect(a:0.1:b)
z = collect(a:0.1:b)


Vboxinfo = VirtualBox(y,z,σ)
N = Vboxinfo.N
dt = 0.01

U₀ = 1.0 #Convective Velocity
#TI #turbulence intensity


#Isotropic turbulence
u_p = (U₀ * TI)^2

Re_stress = [u_p 0.0 0.0; 
            0.0 u_p 0.0;
            0.0 0.0 u_p]



A = cholesky_decomposition(Re_stress)
Eddies = initialize_eddies(N, σ, Vboxinfo)

# Or you can use the wrapper function
A, Eddies = initialize_eddies(U₀, TI, Vboxinfo)


#Computing the velocity in the middle of the VirtualBox domain
vector_points = [[0.0, b/2, b/2]]

#Defining how many time interval

Nt = 1000
U = zeros(Nt, 3)


for i = 1:1:Nt
    U[i,:] = compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo, A)[1]
end

#The deviation standard should approach the turbulence intensity
s1 = Statistics.std(U[:,1])
s2 = Statistics.std(U[:,2])
s3 = Statistics.std(U[:,3])
return mean([s1,s2,s3])
end




