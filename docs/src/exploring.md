# Exploring the package

### Visualize the centre of the eddies
Notice that the default value of the number of eddies is overwritten to reduce the total number in order to make the visualization clearer

```julia
using SyntheticEddyMethod
using Plots
σ = 0.1
b = 5.0
a = 0.0
x = collect(-0.1:0.1:0.1)
y = collect(a:0.1:b)
z = collect(a:0.1:b)


Vboxinfo = VirtualBox(x, y, z, σ)
Vboxinfo.N = 100
dt = 0.01

U₀ = 1.0 #Convective Velocity
TI = 0.2 #turbulence intensity
Re, Eddies = initialize_eddies(U₀, TI, Vboxinfo)

Plots.scatter()
for i = 1:1:length(Eddies)
Plots.scatter!([Eddies[i].xᵢ[1]], [Eddies[i].xᵢ[2]], 
[Eddies[i].xᵢ[3]],   legend=false,  ms=2, color=:black)

end
Plots.scatter!(xlabel="x",ylabel="y",zlabel="z")
```



### Visualize isosurface velocity

```julia
using PlotlyJS

dt = 0.01
X, Y, Z = mgrid(x, y, z)
vector_points = create_vector_points(x, y, z)
Vboxinfo = VirtualBox(x, y, z, σ)
Re, Eddies = initialize_eddies(U₀, TI, Vboxinfo)

Vboxinfo.N

value = map(x-> compute_uSEM(x, Eddies, Vboxinfo, Re)[1], vector_points)

A = value[1]'
for i = 1:1:length(value)
A = vcat(A,value[i]')
end

iso_surfaces = isosurface(
    x=X[:],
    y=Y[:],
    z=Z[:],
    value=A[:,1],
    isomin=0.1,
    isomax=1,
    surface_count=3,
    opacity=0.5,
    caps=attr(x_show=false, y_show=false)
)

layout=Layout(yaxis=attr(scaleanchor="x", scaleratio=1), zaxis=attr(scaleanchor="x", scaleratio=1))
io = PlotlyJS.plot(iso_surfaces, Layout(yaxis=attr(scaleanchor="x", scaleratio=1)))

```

### Visualize the Divergence Free Plane
```julia
using SyntheticEddyMethod
using Statistics
using LinearAlgebra
using Plots
σ = 0.05
b = 1.0
a = 0.0
y = collect(a:0.1:b)
z = collect(a:0.1:b)


Vboxinfo = VirtualBox(y, z, σ)

dt = 0.01
U₀ = 1.0 #Convective Velocity
TI = 0.1 #turbulence intensity
Re, Eddies = initialize_eddies(U₀, TI, Vboxinfo)

#Creating vector of points (0, yp, zp)
val = create_vector_points(0.0,y,z)
dl = 0.0001

#Creating vector of points (dl, yp, zp), (0, yp+dl, zp), (0, yp, zp+dl)
#for computing the divergence 
ux = create_vector_points(dl,y,z)
uy = create_vector_points(0.0, y .+dl,z)
uz = create_vector_points(0.0, y ,z .+ dl)


D = Vector{Float64}[]
for (u,ux,uy,uz) in zip(val, ux,uy,uz)
push!(D,u)
push!(D,ux)
push!(D,uy)
push!(D,uz)
end


r = compute_fluct(D, dt, Eddies, U₀, Vboxinfo::VirtualBox, Re; DFSEM=true)

gD = Float64[]


for i = 1:4:length(D)
dudx = (r[i+1][1] - r[i][1]) / dl

dvdy = (r[i+2][2] - r[i][2]) / dl

dwdz = (r[i+3][3] - r[i][3]) / dl
grad_norm = norm([dudx, dvdy, dwdz])

div_norm =  (dudx + dvdy + dwdz)/norm([dudx, dvdy, dwdz])

if isnan(div_norm)
    push!(gD, 1e-7)
else
    push!(gD, div_norm)
end

end


gD_mat = reshape(gD, (length(z), length(y)))

using LaTeXStrings
tv = -7:-1
tl = [L"10^{%$i}" for i in tv]

p1 = contourf(y, y, log10.(abs.(gD_mat)), fill = true, color=:turbo, aspect_ratio=:equal, levels=8)
plot!(xlabel="y", ylabel="z")
```