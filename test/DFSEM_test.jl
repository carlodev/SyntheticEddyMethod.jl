using SyntheticEddyMethod
using Statistics
using LinearAlgebra
σ = 0.1
b = 5.0
a = 0.0
x = collect(-0.1:0.1:0.1)
y = collect(a:0.1:b)
z = collect(a:0.1:b)


Vboxinfo = VirtualBox(y, z, σ; shape_fun=DFSEM_fun)
dt = 0.01

U₀ = 1.0 #Convective Velocity
TI = 0.1 #turbulence intensity
Re, Eddies = initialize_eddies(U₀, TI, Vboxinfo)
k = compute_kp(Re)

eig_vals = eigvals(Re)


M = I(3).*eig_vals
S = Re * inv(M)

eig_vals, S = eigen(Re)
S * sqrt.(2 .*(k .- eigvals(Re)))

compute_RL(Re)


dt = 0.01



vector_points = create_vector_points(0, b / 2, b / 2)

Nt = 1000
U = zeros(Nt, 3)
r, e = compute_uDFSEM(vector_points, dt, Eddies, U₀, Vboxinfo, Re)

for i = 1:1:Nt
    r, e = compute_uDFSEM(vector_points, dt, Eddies, U₀, Vboxinfo, Re)
    println(i)
    U[i, :] = r[1]
end

U

s1 = Statistics.std(U[:, 1])
s2 = Statistics.std(U[:, 2])
s3 = Statistics.std(U[:, 3])

dl = 0.0001

vv = [[0.0, b / 2, b / 2], [dl, b / 2, b / 2], [0.0, b / 2 + dl, b / 2], [0.0, b / 2, b / 2 + dl]]

r, e = compute_uDFSEM(vv, dt, Eddies, U₀, Vboxinfo, Re)

dudx = (r[2][1] - r[1][1]) / dl

dvdy = (r[3][2] - r[1][2]) / dl

dwdz = (r[4][3] - r[1][3]) / dl

dudx + dvdy + dwdz