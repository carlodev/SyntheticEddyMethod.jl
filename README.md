

# SEM (Synthetic Eddy Method) in Julia

**Documentation**

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://carlodev.github.io/SyntheticEddyMethod.jl/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://carlodev.github.io/SyntheticEddyMethod.jl/)
[![Build Status](https://github.com/carlodev/SyntheticEddyMethod.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/carlodev/SyntheticEddyMethod.jl/actions/workflows/CI.yml?query=branch%3Amaster) 
[![Coverage](https://codecov.io/gh/carlodev/SyntheticEddyMethod.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/carlodev/SyntheticEddyMethod.jl)




## Introduction

The Synthetic Eddy Method (SEM) is a numerical simulation technique used to model turbulent fluid flow in engineering and scientific applications. It involves synthesizing small-scale turbulent structures, or eddies, within a computational domain to represent the effects of larger-scale turbulent flows. This is accomplished by applying perturbations to the flow field, which induce a cascade of energy from larger to smaller eddies until the energy is dissipated through viscous effects. The result is a simulation that captures the important features of turbulent flows while remaining computationally efficient. SEM has been successfully applied to a range of problems, including airfoil and wing design, turbulent combustion, and oceanography. Its ability to accurately capture the physics of turbulent flows makes it a valuable tool for researchers and engineers seeking to improve the efficiency and performance of fluid systems.

The method has been originally developed by Jarrin (10.1016/j.ijheatfluidflow.2006.02.006).

The Divergence-Free Synthetic Eddy Method (DFSEM) is an evolution of the Synthetic Eddy Method (SEM) used for simulating turbulent flows in fluid dynamics. While the SEM uses stochastic generation of eddies to represent the small scales of turbulence, the DFSEM adds the constraint of ensuring that the synthetic eddies produce a divergence-free flow field. In incompressible flows, as the case of turbulent flows, this constraint ensures that the overall flow remains physically consistent and leads to better accuracy and stability in the simulations.

The method has been originally developed by Poletto at al. (10.1007/s10494-013-9488-2).

## Installation
The package is registered, so you can install it as:
```julia
using Pkg
Pkg.add(SyntheticEddyMethod)
```


You can use the most recent release installing it as:
```julia
using Pkg
Pkg.add(url="https://github.com/carlodev/SyntheticEddyMethod.jl")
```

## How to use
At first, the user defines the dimension and resolution of the virtual box where the eddies are going to be generated as well as the dimension of the eddies (`σ`). A common choice is `σ = 2Δz` or `σ = Δz`, where `Δz` is the mesh resolution in the `z` direction.

```julia
using SyntheticEddyMethod

σ = 0.1
b = 5.0
a = 0.0

x = -σ:0.1:+σ
y = collect(a:0.1:b)
z = collect(a:0.1:b)

```

Then you create the `VirtualBox` structure that has embedded all the information about the virtual box where the eddy are generated.
The number of eddy is automatically computed in order to guarantee an homogeneous fill. You can manually override the value (`Vboxinfo.N`).

```julia
Vboxinfo = VirtualBox(y,z,σ)
N = Vboxinfo.N #you can override it 
```

Then, eddies are initialize in the virtualbox with random values of position and intensity. You have to specify the time-step, `dt`. Then the Reynolds stress tensor. Here homegeneous and isotropic turbulence is considered, so $R_{ij} = 0, i!=j; R_{ij} != 0, i=j$, and the terms are computed from the turbulence intensity (`TI`). 

```julia
Eddies = initialize_eddies(Vboxinfo)
t = 0
dt = 0.001
U₀ = 1.0
TI = 0.01 #turbulence intensity

U₀ = 1.0
TI = 0.01 #turbulence intensity

Re_stress, Eddies = initialize_eddies(U₀, TI, Vboxinfo)
```

You have to create a `Vector{Vector{Float64}}` of points where you want to evaluate the speed. In this case you are evaluating the speed in each point of the domain.

```julia
vector_points = create_vector_points(x, y, z)
```
You can create evaluate the speed in just one point (useful for monitoring how the velocity varies in time and creating the spectra)
```julia
vector_points = [[0.0, 1.0, 2.5]]
```

Compute the velocity fluctuation. It is then is 'corrected' using the matrix A which is internally created using the `cholesky_decomposition` function.
```julia
u_fluct = compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo, Re_stress)
```



## Package Features
- Create fluctuations that respect the divergence-free condition (DFSEM)
- Create velocity fluctuations for inlet boundary conditions
- Create coeherent eddies in 3D domain
- Define custom Reynolds Stress Tensor
- Import from file custom Reynolds Stress Tensor

## Acknowledgement
- Nomenclature: 10.1016/j.ijheatfluidflow.2006.02.006
- Shape function definition and Reynolds stress tensor sheet https://github.com/blackcata/SEM.git and the related paper 10.1016/j.ijheatmasstransfer.2019.02.061
- For detailed description of the procedure https://nheri-simcenter.github.io/WE-UQ-Documentation/common/technical_manual/desktop/WEUQ/TinF.html 
- DFSEM: 10.1007/s10494-013-9488-2