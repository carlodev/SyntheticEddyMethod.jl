module SyntheticEddyMethod

using Plots
using DataFrames
using LinearAlgebra
using Statistics
using FFTW #For fast fourier transformation

export tent_fun
export step_fun
export trunc_gauss_fun
export fσ
include("Shapefunctions.jl")

export cholesky_decomposition
include("Decomposition.jl")

export SEM_EDDY
export VirtualBox
export convect_eddy
export create_vector_points
export compute_U_k
export initialize_eddies
export compute_uᵢₚ
include("Eddies.jl")

export fft_from_signal
include("SpectralUtilities.jl")
end
