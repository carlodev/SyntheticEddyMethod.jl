using SyntheticEddyMethod
using Statistics
using LinearAlgebra
using Plots
σ = 0.1
b = 5.0
a = 0.0
x = collect(-0.1:0.1:0.1)
y = collect(a:0.1:b)
z = collect(a:0.1:b)


Vboxinfo = VirtualBox(y, z, σ;)
dt = 0.01

U₀ = 1.0 #Convective Velocity
TI = 0.1 #turbulence intensity
Re, Eddies = initialize_eddies(U₀, TI, Vboxinfo)
vector_points = create_vector_points(0, b / 2, b / 2)
dt = 0.01

uip, e = compute_uSEM(vector_points, dt, Eddies, U₀, Vboxinfo, Re)

uip

Nt = 1000
U = zeros(Nt, 3)
r, e = compute_uSEM(vector_points, dt, Eddies, U₀, Vboxinfo, Re)
r



for i = 1:1:Nt
    r = compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo, Re)
    println(i)
    U[i, :] = r[1]
end

U

s1 = Statistics.std(U[:, 1])
s2 = Statistics.std(U[:, 2])
s3 = Statistics.std(U[:, 3])

