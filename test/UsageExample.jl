using SyntheticEddyMethod

using Statistics
using DataFrames, XLSX
#using PlotlyJS
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
TI = 0.1 #turbulence intensity


#Isotropic turbulence
u_p = (U₀ * TI)^2

Re_stress = [u_p 0.0 0.0; 
            0.0 u_p 0.0;
            0.0 0.0 u_p]



A = cholesky_decomposition(Re_stress)
Eddies = initialize_eddies(Vboxinfo)

# Or you can use the wrapper function
A, Eddies = initialize_eddies(U₀, TI, Vboxinfo)


#Computing the velocity in the middle of the VirtualBox domain
vector_points = [[0.0, b/2, b/2]]

#Defining how many time interval
Nt = 200
q = zeros(Nt, 3)

for i = 1:1:Nt
    q[i,:] = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]
end

U, Ek =  compute_U_k(q, A, U₀)


Statistics.std(U)
Statistics.std(Ek)





## Plotting 3D iso curves (good for visualizing the distribution and evolution of the eddies)

# X, Y, Z = mgrid(x, y, z)
# vector_points = create_vector_points(x, y, z)

# value = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]
# Eddies

# iso_surfaces = isosurface(
#     x=X[:],
#     y=Y[:],
#     z=Z[:],
#     value=value[:,1],
#     isomin=0.1,
#     isomax=1,
#     surface_count=3,
#     opacity=0.5,
#     caps=attr(x_show=false, y_show=false)
# )

# layout=Layout(yaxis=attr(scaleanchor="x", scaleratio=1), zaxis=attr(scaleanchor="x", scaleratio=1))
# io = PlotlyJS.plot(iso_surfaces, Layout(yaxis=attr(scaleanchor="x", scaleratio=1)))