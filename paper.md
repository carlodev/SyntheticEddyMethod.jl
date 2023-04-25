---
title: 'SyntheticEddyMethod.jl: A Julia package for the creation of inlet flow conditions for LES'
tags:
  - Julia
  - Turbulence
  - Eddy
  - Inlet Condtions
  - LES
authors:
  - name: Carlo Brunelli
    orcid: 0000-0002-2873-6293
    equal-contrib: true
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
  - name: Bart Janssens
    corresponding: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Georg May
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3
  - name:  Mark Runacres
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3    
affiliations:
 - name: Mechanical Engineering Department, Royal Military Academy, Belgium
   index: 1
 - name: Aeronautics and Aerospace Department, von Karman Institute for Fluid Dynamics, Belgium
   index: 2
 - name: Engineering Technology Thermodynamics and Fluid Mechanics Group, VUB, Belgium
   index: 3
date: 25 April 2023
bibliography: paper.bib

---

# Summary
The Synthetic Eddy Method (SEM) is a numerical simulation technique used to model turbulent fluid flow in engineering and scientific applications. It involves synthesizing small-scale turbulent structures, or eddies, within a computational domain to represent the effects of larger-scale turbulent flows. This is accomplished by applying perturbations to the flow field, which induce a cascade of energy from larger to smaller eddies until the energy is dissipated through viscous effects. The result is a simulation that captures the important features of turbulent flows while remaining computationally efficient. Its ability to accurately capture the physics of turbulent flows makes it a valuable tool for researchers and engineers seeking to improve the efficiency and performance of fluid systems.

# Statement of need

`SyntheticEddyMethod.jl` is a package which aims to create realistic turbulent inlet conditions for Large Eddies Simulations. The method has been originally introduced by `@Jarrin:2006`. It has been extended by `@Poletto:2013` implementing the divergence-free constraint to fluctuations for incompressible flows.

`SyntheticEddyMethod.jl` is completly implemented in Julia programming language, `@Bezanson:2017`. In recent years, Julia has emerged as a powerful language for scientific computing and has become popular among researchers and practitioners in the field of fluidynamics. Julia is is extremely expressive and allows to condensate complex mathematical expression in few synthetic lines. This package will be a valuable tool for researchers and engineers working in the field of CFD, offering an intuitive and efficient way to simulate proper boundaries conditions.



# Usage Example

```julia
using SyntheticEddyMethod
σ = 0.1 #eddy dimension

#VirtualBox size
b = 5.0
a = 0.0

x = -σ:0.1:+σ
y = collect(a:0.1:b)
z = collect(a:0.1:b)

#VirtualBox creation
Vboxinfo = VirtualBox(y,z,σ)

dt = 0.001 #time step
U₀ = 1.0 #convective velocity
TI = 0.01 #turbulence intensity

#Creating eddies and Reynolds stress tensor
Re_stress, Eddies = initialize_eddies(U₀, TI, Vboxinfo) 

eval_point = [0.0, 1.0, 2.5] #evaluation point

u_fluct = compute_fluct(eval_point, dt, Eddies, U₀, Vboxinfo, Re_stress)

```

This is an example of the spectra created using SyntheticEddyMethod with a tent function for different turbulence intensity. The spectra referers to the fluctutation in time in one speicific point.

![Spectra](images/docs/Spectra.png){ width=50% }


It is reported the normalized divergence in a plane using the DFSEM.

![Divergence Free](images/docs/Div_free_plane.png){ width=50% }

# Package Features
- Create fluctuations that respect the divergence-free condition (DFSEM)
- Create velocity fluctuations for inlet boundary conditions
- Create coeherent eddies in 3D domain
- Define custom Reynolds Stress Tensor
- Import from file custom Reynolds Stress Tensor


# References
