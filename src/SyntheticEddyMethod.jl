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
export DFSEM_fun
export fσ
include("Shapefunctions.jl")

export cholesky_decomposition
include("Decomposition.jl")

export Reynolds_stress_interpolator
export Reynolds_stress_tensor
export get_reynolds_stress_from_file
include("ReynoldsStress.jl")

export AbstractEddy
export SemEddy
export VirtualBox
export convect_eddy
export initialize_eddies
include("Eddies.jl")

export compute_fluct
export compute_uSEM
export compute_uDFSEM
export compute_rk
export compute_RL
export compute_α
export compute_kp
include("Fluctuations.jl")

export create_vector_points
export compute_Ek
export fft_from_signal
include("Utilities.jl")

end

