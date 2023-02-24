# Package usage
At first, the user defines the dimension and resolution of the virtual box where the eddies are going to be generated as well as the dimension of the eddies (`σ`). A common choice is `σ = 2Δz` or `σ = Δz`, where `Δz` is the mesh resolution in the `z` direction.

```julia
using SyntetichEddyMethod

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
You can specify also the x direction dimensions, a different shape funcion and different ´σ´ in the three different directions.

```julia
using SyntetichEddyMethod

σ = 0.1
b = 5.0
a = 0.0
a = collect(-1:0.1:1)
y = collect(a:0.1:b)
z = collect(a:0.1:b)
σx = 0.1
σy = 0.05
σz = 0.07
σ = [σx,σy,σz]
Vboxinfo = VirtualBox(x,y,z,σ; shape_fun = step_fun)
```



Then, eddies are initialize in the virtualbox with random values of position and intensity. You have to specify the time-step, `dt`. Then the Reynolds stress tensor. Here homegeneous and isotropic turbulence is considered, so $R_{ij} = 0, i!=j; R_{ij} != 0, i=j$, and the terms are computed from the turbulence intensity (`TI`). Then the Matrix `A` is created using the `cholesky_decomposition` function.


```julia
Eddies = initialize_eddies(N, σ, Vboxinfo)
t = 0
dt = 0.001
U₀ = 1.0
TI = 0.6 #turbulence intensity

#Isotropic turbulence
u_p = (U₀ * TI)^2

Re_stress = [u_p 0.0 0.0; 
            0.0 u_p 0.0;
            0.0 0.0 u_p]

A = cholesky_decomposition(Re_stress)

```

You have to create a `Vector{Vector{Float64}}` of points where you want to evaluate the speed.

```julia
vector_points = create_vector_points(x, y, z)
```
You can create evaluate the speed in just one point (useful for monitoring how the velocity varies in time and creating the spectra)
```julia
vector_points = [[0.0, 1.0, 2.5]]
```

Compute the velocity fluctuation and then is 'corrected' using the matrix A. The velocity and the turbulent kinetic energy are then computed.

```julia
u_fluct = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]
U, Ek =  compute_U_k(u_fluct, A, U₀)
```

From here you can generate the spectra of your signal. For more detail look at files in test/


## Tests
Look in folder /test for the in-deep results of testing




