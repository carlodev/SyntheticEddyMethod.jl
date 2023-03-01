module SyntheticEddyMethod

using Plots
using DataFrames
using LinearAlgebra
using Statistics
using FFTW #For fast fourier transformation
using Interpolations
using DataFrames
using XLSX


export tent_fun
export step_fun
export trunc_gauss_fun
export fσ
include("Shapefunctions.jl")

export cholesky_decomposition
include("Decomposition.jl")

export Reynolds_stress_interpolator
export Reynolds_stress_tensor
export Reynolds_stress_points
export Reynolds_stress_point
export get_reynolds_stress_from_file
include("ReynoldsStress.jl")

export AbstractEddy
export SemEddy
export VirtualBox
export convect_eddy
export create_vector_points
export initialize_eddies
export compute_fluct
export compute_uᵢₚ
export compute_Ek
export uᵢ
export new_rand_position
include("Eddies.jl")



export fft_from_signal
include("SpectralUtilities.jl")
end
