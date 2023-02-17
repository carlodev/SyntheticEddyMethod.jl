
using Statistics
function driver_test(TI::Float64; Nt = 20000)
σ = 0.1 #eddy dimensions, the same in all the directions
b = 5.0
a = 0.0

#Defining the Virtual Box domain
x = -σ:0.1:+σ 
y = a:0.1:b
z = a:0.1:b


Vboxinfo = VirtualBox(y,z,σ)

N = Vboxinfo.N #you can override it 
t = 0
dt = 0.001

U₀ = 1.0 #Convective Velocity
#TI = 0.1 #turbulence intensity


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
#Nt = 20000
q = zeros(Nt, 3)

for i = 1:1:Nt
    q[i,:] = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]
end

U, Ek =  compute_U_k(q, A, U₀)

return Statistics.std(U)

end




