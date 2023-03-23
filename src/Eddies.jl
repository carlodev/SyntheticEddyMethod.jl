"""
    AbstractEddy

Abstract type defining an Eddy object.
"""
abstract type AbstractEddy end

"""
    SemEddy

It identify the properties of each eddy.
 - eddy_num::Int64     : Eddy identification number
 - σ::Vector{Float64}  : Eddy length scale
 - xᵢ::Vector{Float64} : Eddy's position in the computational box [x,y,z]
 - ϵᵢ::Vector{Float64}  : Eddy's intensity (+1 or -1) in [x,y,z]
"""
mutable struct SemEddy <: AbstractEddy
    eddy_num::Int64     # Eddy identification number
    σ::Vector{Float64}  # Eddy length scale
    xᵢ::Vector{Float64} # Eddy's position in the computational box [x,y,z]
    ϵᵢ::Vector{Float64}  # Eddy's intensity (+1 or -1) in [x,y,z]
end


"""
    VirtualBox

Virtual Volume box where the eddies are created.
 - σ::Vector{Float64} : 3 elements vector, to have different σ in different directions
 - N::Int64 : number of eddies
 - shape_fun::Function : shape function
 - V_b::Float64 : volume of the virtual box
 - X_start::Float64
 - X_end::Float64
 - Y_start::Float64
 - Y_end::Float64
 - Z_start::Float64
 - Z_end::Float64
"""
mutable struct VirtualBox
    σ::Vector{Float64} # 3 elements vector, to have different σ in different directions
    N::Int64
    T::Float64
    shape_fun::Function
    V_b::Float64
    X_start::Float64
    X_end::Float64
    Y_start::Float64
    Y_end::Float64
    Z_start::Float64
    Z_end::Float64
end

#Virtual Box constructor
"""
If not specified, in the x direction the dimension in from -σ to + σ
"""
function VirtualBox(Y::Vector{Float64}, Z::Vector{Float64},  σ::Float64; shape_fun = tent_fun)
    X = [0.0]
    VirtualBox(X, Y, Z, σ; shape_fun = shape_fun)
end

function VirtualBox(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64},  σ::Float64; shape_fun = tent_fun)
    VirtualBox(X, Y, Z, [σ, σ, σ]; shape_fun = shape_fun)
end

function VirtualBox(X::Vector{Float64}, Y::Vector{Float64}, Z::Vector{Float64},  σ::Vector{Float64}; shape_fun = tent_fun)
    X_start = X[1] - σ[1]
    X_end = X[end] + σ[1]
    Y_start = Y[1] - σ[2]
    Y_end = Y[end] + σ[2]
    Z_start = Z[1] - σ[3]
    Z_end = Z[end] + σ[3]
    #Define a suitable Law For defining Sp for longer domains
    Sₚ = (Y_end - Y_start) * (Z_end - Z_start) * ((X_end - X_start)/ (2*σ[1]))^(1)
    Sₛ = σ[2] * σ[3]  #Eddy surface on the XY Plane
    N = Int(round(Sₚ / Sₛ))
    V_b =  (Y_end - Y_start) * (Z_end - Z_start) * (X_end - X_start)
    T = 0.0 #Set the initialization time
    VirtualBox(σ, N, T, shape_fun, V_b, X_start, X_end, Y_start, Y_end, Z_start, Z_end)
end



"""
    initialize_eddies(Vbinfo::VirtualBox)

Initialize Eddy position and intensity
"""
function initialize_eddies(Vbinfo::VirtualBox)
    SEMEddy = SemEddy[]
    for i =1:1:Vbinfo.N
        ϵᵢ = rand((-1,1), 3)
        xᵢ = new_rand_position(Vbinfo)
        push!(SEMEddy, SemEddy(i, Vbinfo.σ,  xᵢ,  ϵᵢ))
    end
    return SEMEddy
end


function initialize_eddies(U₀::Real, TI::Float64, Vboxinfo::VirtualBox; turbulence_type= :hom_is)

    u_p = (U₀ * TI)^2
    
    if turbulence_type == :hom_is #homogeneous and isotropic
        Re_stress = collect(I(3).*u_p)
    
    else
        "Here if you want to implement different type of turbulence"
        error("Turbulence type supported :hom_is")
    end
    
    Eddies = initialize_eddies(Vboxinfo)
    
    return Re_stress, Eddies    
end


# 
#     new_rand_position(Vbinfo::VirtualBox)

# It computes a random position inside the Virtual Box volume.
# 
function new_rand_position(Vbinfo::VirtualBox)

    xx = (rand() .- 0.5) .* (Vbinfo.X_end - Vbinfo.X_start) .+ (Vbinfo.X_end + Vbinfo.X_start) ./ 2 #Vbinfo.X_start 
    yy = (rand() .- 0.5) .* (Vbinfo.Y_end - Vbinfo.Y_start) .+ (Vbinfo.Y_end + Vbinfo.Y_start) ./ 2
    zz =(rand() .- 0.5) .* (Vbinfo.Z_end - Vbinfo.Z_start) .+ (Vbinfo.Z_end + Vbinfo.Z_start) ./ 2

    return[xx,yy,zz]
end


"""
    convect_eddy(dt::Float64, Eddy::SemEddy, U₀::Float64, σ::Vector{Float64}, Vbinfo::VirtualBox)

The eddies are convected, shifted in the `x` direction of a distance `dt * U₀`.\\
If one eddy exits the Virtual Box, it is re-generated randomly inside the Virtual Box.
"""
function convect_eddy(dt::Float64, Eddy::SemEddy, U₀::Float64, Vbinfo::VirtualBox)
    x_tmp = Eddy.xᵢ[1] + dt * U₀
    if x_tmp < Vbinfo.X_end
        Eddy.xᵢ = [x_tmp, Eddy.xᵢ[2], Eddy.xᵢ[3]]
    else
        Eddy.xᵢ = new_rand_position(Vbinfo)
        Eddy.ϵᵢ = rand((-1,1), 3)
    end
       
    return Eddy
end