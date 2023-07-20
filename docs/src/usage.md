# Package usage

## Using the SEM
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
Then, eddies are initialize in the virtualbox with random values of position and intensity. You have to specify the time-step, `dt`. The Reynolds stress tensor can be specified by the user (as a 3x3 matrix) or automatically computed just providing the turbulence intensity and convective velocity.


```julia
Eddies = initialize_eddies(Vboxinfo)
t = 0
dt = 0.001
U₀ = 1.0 #Convective velocity, x-axis
TI = 0.01 #turbulence intensity

Re_stress, Eddies = initialize_eddies(U₀, TI, Vboxinfo)

```

You have to create a `Vector{Vector{Float64}}` of points where you want to evaluate the speed.
```julia
x = 0.0
vector_points = create_vector_points(x, y, z)

```

```julia
u_fluct = compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo, Re_stress)
```

You can create evaluate the speed in just one point (useful for monitoring how the velocity varies in time and creating the spectra)
```julia
single_point = [0.0, 1.0, 2.5]
u_fluct = compute_fluct(single_point, dt, Eddies, U₀, Vboxinfo, Re_stress)

```

Compute the velocity fluctuation and then is 'corrected' using the Reynolds Stress tensor.



Compute the turbulent kinetic energy:
```julia
    Ek  = compute_Ek(u_fluct, U₀)
```

## Using the DFSEM
Following an analogus procedure is possible to use the divergence-free sythetic eddy method. It allows to create fluctuations that are divergence-free, useful for incompressible flows.
The virtual box is created in an analogus way. Internally the shape function in this case is [`DFSEM_fun`](@ref) which has been specifically designed for the DFSEM.
```julia
using SyntheticEddyMethod
using LinearAlgebra
σ = 0.1
b = 5.0
a = 0.0
y = collect(a:0.1:b)
z = collect(a:0.1:b)

Vboxinfo = VirtualBox(y,z,σ)
dt = 0.001
U₀ = 1.0
TI = 0.01 #turbulence intensity

Re_stress, Eddies = initialize_eddies(U₀, TI, Vboxinfo)

```
Now, we are going to evaluate the fluctuation in 4 points useful for approximating the divergence:\
``[x,y,z], [x+dx, y, z], [x,y+dy,z], [x,y,z+dz]``

```julia
dl = 0.0001
vector_points = [[0.0, b / 2, b / 2], [dl, b / 2, b / 2], [0.0, b / 2 + dl, b / 2], [0.0, b / 2, b / 2 + dl]]
u_fluct = compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo, Re_stress; DFSEM = true)
```

In order to verify the diverge-free condition is respected, the derivatives in each direction are computed with a simple forward method.
The divergence is normalized with the module of the gradient to obtain a non-diemensional quantity.\
``\dfrac{\nabla\cdot \vec{u}}{|\nabla \vec{u}|}``
```julia
dudx = (u_fluct[2][1] - u_fluct[1][1]) / dl
dvdy = (u_fluct[3][2] - u_fluct[1][2]) / dl
dwdz = (u_fluct[4][3] - u_fluct[1][3]) / dl
grad_norm = norm([dudx, dvdy, dwdz])
div_val = dudx + dvdy + dwdz

div_val/grad_norm
```
Note: if ``|\nabla \vec{u}|`` is too small, it can leads `NaN` values.




