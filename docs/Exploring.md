# Exploring the package

### Visualize the centre of the eddies
```julia
using SyntheticEddyMethod
using Plots
σ = 0.1
b = 5.0
a = 0.0
x = collect(-0.1:0.1:0.1)
y = collect(a:0.1:b)
z = collect(a:0.1:b)


Vboxinfo = VirtualBox(x, y,z, σ)
Vboxinfo.N = 100
dt = 0.01

U₀ = 1.0 #Convective Velocity
TI = 0.2 #turbulence intensity
A, Eddies = initialize_eddies(U₀, TI, Vboxinfo)

Plots.scatter([Eddies[1].xᵢ[1]], [Eddies[1].xᵢ[2]], [Eddies[1].xᵢ[3]], legend=false, ms=2, color=:black)
for i = 2:1:length(Eddies)-1
    Plots.scatter!([Eddies[i].xᵢ[1]], [Eddies[i].xᵢ[2]] , [Eddies[i].xᵢ[3]],  legend=false,  ms=2, color=:black)

end
Plots.scatter!([Eddies[end].xᵢ[1]], [Eddies[end].xᵢ[2]] , [Eddies[end].xᵢ[3]],   legend=false,  ms=2, color=:black)


```



### Visualize isosurface velocity
```julia
using PlotlyJS
X, Y, Z = mgrid(x, y, z)
vector_points = create_vector_points(x, y, z)

value = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]

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
```

