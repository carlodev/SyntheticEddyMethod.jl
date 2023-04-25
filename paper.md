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
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: Mechanical Engineering Department, Royal Military Academy, Belgium
   index: 1
date: 25 April 2023
bibliography: paper.bib

---

# Summary
The Synthetic Eddy Method (SEM) is a numerical simulation technique used to create turbulent flow field with desired features. It is used in computational fluid dynamics for imposing realistic inlet boundary conditions improving the fidelity of the results obtained by simulations. Its ability to create fluctuation with prescribed physical features makes it a valuable tool for researchers and engineers seeking to improve the reliability of simulations and also for reacreating as much as possible an evironment close to an experimental one. The package allows users to easily generate synthetic turbulence fields that can be used in CFD simulations, and to control the level of turbulence and eddies length of the generated fields.

# Statement of need

`SyntheticEddyMethod.jl` is a package which aims to create realistic turbulent inlet conditions for Large Eddies Simulations. This package will be a valuable tool for researchers and engineers working in the field of Computational Fluid Dynamics, offering an intuitive and efficient way to simulate proper boundary conditions. The fluctuations generated are more realistic than those that can be easily produced by a random signal.

The method has been originally introduced by @Jarrin:2006. It has been extended by @Poletto:2013 implementing the divergence-free (DFSEM) constraint to fluctuations for incompressible flows. 

`SyntheticEddyMethod.jl` is completly implemented in Julia programming language, [@Bezanson:2017]. In recent years, Julia has emerged as a powerful language for scientific computing and has become popular among researchers and practitioners in the field of fluid dynamics. Julia is extremely expressive and allows to condensate complex mathematical expression in a few synthetic lines. The functions are written almost identically as on paper. This has also an advantage for people who desire to contribute and use it. By taking advantage of the flexibility of Julia multiple dispatch, it allows users to simulate fluctuations at specific points in the flow field or at multiple points simultaneously, offering a powerful optimized tool.Users can customize several key parameters of the SEM method, such as the turbulence intensity, Reynolds stress, and eddy dimensions. These parameters can be set by the user directly, or loaded from a file, making the package versatile and user-friendly. Additionally, it can simulate anisotropic effects by allowing the eddies to have different dimensions along different directions, enabling users to model complex flows accurately. 

Julia allows the package to be flexible, for example computing the fluctuations in one point of the domain or in multiple locations or providing one Reynolds stress tensor the whole domain or pointwise. 

Different software packages have been developed to implement this method (for example using Fortran @OH:2019). However, these packages are often limited in their applicability and can be challenging for non-experts to use. `SyntheticEddyMethod.jl` is designed to be more general-purpose, allowing it to be applied to a broader range of turbulence simulation problems. It is designed to be more accessible and with clear documentation.


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

## Results
![Spectra. \label{Spectra}](images/docs/Spectra.png){ width=50% }
![DFSEM plane. \label{dfsemPlane}](images/docs/Div_free_plane.png){ width=50% }

Figure \autoref{Spectra} Example of the spectra created using SyntheticEddyMethod with a tent function for different turbulence intensities. The spectra in figure \autoref{dfsemPlane} refers to the fluctuations in time in one specific point. Normalized divergence in a plane using the divergence-free Synthetic Eddy Method.


# Package Features
- Create fluctuations that respect the divergence-free condition (DFSEM)
- Create velocity fluctuations for inlet boundary conditions
- Create coeherent eddies in 3D domain
- Define custom Reynolds Stress Tensor
- Import from file custom Reynolds Stress Tensor


# References