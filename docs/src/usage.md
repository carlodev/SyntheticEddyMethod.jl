# Package usage
At first, the user defines the dimension and resolution of the virtual box where the eddies are going to be generated as well as the dimension of the eddies (`σ`). A common choice is `σ = 2Δz` or `σ = Δz`, where `Δz` is the mesh resolution in the `z` direction.

```julia
using SyntheticEddyMethod

σ = 0.1
b = 5.0
a = 0.0

y = collect(a:0.1:b)
z = collect(a:0.1:b)

```
Notice: only the first and last value of ´y´ and ´z´ are used to create the VirtualBox.

Then you create the `VirtualBox` structure that has embedded all the information about the virtual box where the eddy are generated.
The number of eddy is automatically computed in order to guarantee an homogeneous fill. You can manually override the value (`Vboxinfo.N`).

```julia

Vboxinfo = VirtualBox(y,z,σ)

N = Vboxinfo.N #you can override it 
```
Then, eddies are initialize in the virtualbox with random values of position and intensity. You have to specify the time-step, `dt`. Then the Reynolds stress tensor. Here homegeneous and isotropic turbulence is considered, so $R_{ij} = 0, i!=j; R_{ij} != 0, i=j$, and the terms are computed from the turbulence intensity (`TI`). Then the Matrix `A` is created using the `cholesky_decomposition` function.


```julia
Eddies = initialize_eddies(N, σ, Vboxinfo)
t = 0
dt = 0.001
U₀ = 1.0
TI = 0.01 #turbulence intensity

A, Eddies = initialize_eddies(U₀, TI, Vboxinfo)

```

You have to create a `Vector{Vector{Float64}}` of points where you want to evaluate the speed.
```julia
x = 0.0
vector_points = create_vector_points(x, y, z)

```
You can create evaluate the speed in just one point (useful for monitoring how the velocity varies in time and creating the spectra)
```julia
vector_points = [[0.0, 1.0, 2.5]]
```

Compute the velocity fluctuation and then is 'corrected' using the Reynolds Stress tensor.


```julia
u_fluct = compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo,A )
```

Compute the turbulent kinetic energy
```julia
    Ek  = compute_Ek(u_fluct, U₀)
```


From here you can generate the spectra of your signal. For more detail look at files in test/





