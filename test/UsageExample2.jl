using Revise

using SyntheticEddyMethod

using Statistics
using DataFrames, XLSX
σ = 0.1 #eddy dimensions, the same in all the directions
b = 5.0
a = 0.0

#Defining the Virtual Box domain
x = collect(-σ:0.1:σ) 
y = collect(a:0.1:b)
z = collect(a:0.1:b)

Vboxinfo = VirtualBox(y,z,σ;shape_fun = tent_fun)
Vboxinfo.V_b


N = Vboxinfo.N #you can override it 
t = 0
dt = 0.01
U₀ = 1.0 #Convective Velocity
TI = 0.1 #turbulence intensity


#Isotropic turbulence
u_p = (U₀ * TI)^2

Re_stress = [u_p 0.0 0.0; 
            0.0 u_p 0.0;
            0.0 0.0 u_p]



A = cholesky_decomposition(Re_stress)
Eddies = initialize_eddies(N, [σ,σ,σ], Vboxinfo)

# Or you can use the wrapper function
A, Eddies = initialize_eddies(U₀, TI, Vboxinfo)


#Computing the velocity in the middle of the VirtualBox domain
vector_points = [[0.0, b/2, b/2]]

# Defining how many time interval
 NN = 100
 Stat = zeros(NN,3)

for j = 1:1:NN
    println(j)
Nt = 100_000
q = zeros(Nt, 3)

for i = 1:1:Nt
    q[i,:] = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]
end

U, Ek =  compute_U_k(q, A, U₀)


Stat[j,1] = Statistics.std(U[:,1])
Stat[j,2] = Statistics.std(U[:,2])
Stat[j,3] = Statistics.std(U[:,3])
end

using FileIO, JLD2
Stat[1:50,:]
FileIO.save("StatF2.jld2","Stat",Stat[1:50,:])
FileIO.save("UF2.jld2","U",U)


# Plotting 3D iso curves (good for visualizing the distribution and evolution of the eddies)

using Plotly, Plots
σ = 0.1 #eddy dimensions, the same in all the directions
b = 5.0
a = 0.0


x = collect(-2*σ:0.02:2*σ) 
y = collect(a:0.05:b)
z = [0.0]


Vboxinfo = VirtualBox(x, y,z,σ;shape_fun= step_fun)


N = Vboxinfo.N #you can override it 
t = 0
dt = 0.01

U₀ = 1.0 #Convective Velocity
TI = 0.2 #turbulence intensity
A, Eddies = initialize_eddies(U₀, TI, Vboxinfo)

Plots.scatter([Eddies[1].xᵢ[1]], [Eddies[1].xᵢ[2]],  legend=false, ms=2, color=:black)
for i = 2:1:length(Eddies)-1
    Plots.scatter!([Eddies[i].xᵢ[1]], [Eddies[i].xᵢ[2]],  legend=false,  ms=2, color=:black)

end
Plots.scatter!([Eddies[end].xᵢ[1]], [Eddies[end].xᵢ[2]],  legend=false,  ms=2, color=:black)


vector_points = create_vector_points(x, y, z)


value = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]
vv = value[:,1]

Z = zeros(length(y), length(x))
for i = 1:1:length(y)
    start_idx = length(x)*(i-1)+1
    end_idx = length(x)*i
    Z[i,:]= vv[start_idx:end_idx]
end


X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(y, 1, length(x))
plotly()
p1 = Plots.contour(x, y, Z, fill = true, aspect_ratio=:equal)

p2 = Plots.contour(x, y, Z, aspect_ratio=:equal)
plot(p1, p2)

Plots.savefig(p1,"eddy.html")

X, Y, Z = mgrid(x, y, z)
vector_points = create_vector_points(x, y, z)

value = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]


plotlyjs()
iso_surfaces = isosurface(
    x=X[:],
    y=Y[:],
    z=Z[:],
    value=value[:,1],
    isomin=0.1,
    isomax=1,
    surface_count=3,
    opacity=0.5,
    caps=attr(x_show=false, y_show=false)
)

layout=Layout(yaxis=attr(scaleanchor="x", scaleratio=1), zaxis=attr(scaleanchor="x", scaleratio=1))
io = PlotlyJS.plot(iso_surfaces, Layout(yaxis=attr(scaleanchor="x", scaleratio=1)))



x = collect(1:5)
y = collect(1:2)
z = [0.0]
sigma_vec = [0.1,0.2,0.3]
ck = [1.0,2.0,3.0]
vector_points = create_vector_points(x, y, z)


s2 = map(y -> y ./sigma_vec, map(x -> x .- ck, vector_points))
