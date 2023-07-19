# SyntheticEddyMethod.jl
*Documentation of SyntheticEddyMethod.jl for synthetic eddy generation*

## Introduction
The Synthetic Eddy Method (SEM) is a numerical simulation technique used to model turbulent fluid flow in engineering and scientific applications. It involves synthesizing small-scale turbulent structures, or eddies, within a computational domain to represent the effects of larger-scale turbulent flows. This is accomplished by applying perturbations to the flow field, which induce a cascade of energy from larger to smaller eddies until the energy is dissipated through viscous effects. The result is a simulation that captures the important features of turbulent flows while remaining computationally efficient. SEM has been successfully applied to a range of problems, including airfoil and wing design, turbulent combustion, and oceanography. Its ability to accurately capture the physics of turbulent flows makes it a valuable tool for researchers and engineers seeking to improve the efficiency and performance of fluid systems.

The method has been originally developed by [Jarrin2006](@cite).

The Divergence-Free Synthetic Eddy Method (DFSEM),originally developed by [Poletto2013](@cite), is an evolution of the Synthetic Eddy Method (SEM) used for simulating turbulent flows in fluid dynamics. While the SEM uses stochastic generation of eddies to represent the small scales of turbulence, the DFSEM adds the constraint of ensuring that the synthetic eddies produce a divergence-free flow field. In incompressible flows, as the case of turbulent flows, this constraint ensures that the overall flow remains physically consistent and leads to better accuracy and stability in the simulations.


## Package Features
- Create fluctuations that respect the divergence-free condition (DFSEM)
- Create velocity fluctuations for inlet boundary conditions
- Create coherenteddies in 3D domain
- Define custom Reynolds Stress Tensor
- Import from file custom Reynolds Stress Tensor

## Acknowledgement
- nomenclature: [Jarrin2006](@cite)
- shape function definition thanks to the Fortran 90 code https://github.com/blackcata/SEM.git and the related paper [OH2019](@cite)
- https://nheri-simcenter.github.io/WE-UQ-Documentation/common/technical_manual/desktop/WEUQ/TinF.html for detailed description of the procedure