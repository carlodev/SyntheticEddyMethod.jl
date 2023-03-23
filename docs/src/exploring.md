# Exploring the package

### Visualize the centre of the eddies
Notice that the default value of the number of eddies is overwritten to reduce the total number and make the visualization easier

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

