"""
Structure with the property of a single eddy
"""
mutable struct SEM_EDDY
    eddy_num::Int64     # Eddy specification number
    σ::Float64  # Eddy length scale
    xᵢ::Vector{Float64} # Eddy's position in the computational box [x,y,z]
    ϵᵢ::Vector{Float64}  # Eddy's intensity (+1 or -1) in [x,y,z]

end


"""
Virtual Volume box where the eddies are created
"""
struct VirtualBox
    Y::Vector{Float64}
    Z::Vector{Float64}
    σ::Float64

    N::Int64
    V_b::Float64
    Y_start::Float64
    Y_end::Float64
    Z_start::Float64
    Z_end::Float64

    function VirtualBox(Y, Z, σ)
        Y_start = Y[1] - σ
        Y_end = Y[end] + σ
        Z_start = Z[1] - σ
        Z_end = Z[end] + σ
  

        Sₚ = (Y_end - Y_start) * (Z_end - Z_start)
        Sₛ = σ * σ #Eddy surface on the XY Plane
        N = Int(round(Sₚ / Sₛ))
        V_b = 2*σ * (Y_end - Y_start) * (Z_end - Z_start)

        new(Y, Z, σ, N, V_b, Y_start, Y_end, Z_start, Z_end)
    end

end



"""
Initialize Eddy position and intensity
"""
function initialize_eddies(N::Int64, σ::Float64, Vbinfo::VirtualBox)
    SEM_Eddy = SEM_EDDY[]
    for i =1:1:N
        ϵᵢ = rand((-1,1), 3)
        xᵢ = new_rand_position(Vbinfo)
        push!(SEM_Eddy, SEM_EDDY(i, σ,  xᵢ,  ϵᵢ))
    end
    return SEM_Eddy
end


function initialize_eddies(U₀::Real, TI::Float64, Vboxinfo::VirtualBox; turbulence_type= :hom_is)

    u_p = (U₀ * TI)^2
    
    if turbulence_type == :hom_is #homogeneous and isotropic
    Re_stress = [u_p 0.0 0.0; 
                0.0 u_p 0.0;
                0.0 0.0 u_p]
    
    else
        "Here if you want to implement different type of turbulence"
        error("Turbulence type supported :hom_is")
    end
    
    A = cholesky_decomposition(Re_stress)
    Eddies = initialize_eddies(Vboxinfo.N, Vboxinfo.σ, Vboxinfo)
    
    return A, Eddies    
end

"""
Random position of an eddy inside the Virtual Box
"""
function new_rand_position(Vbinfo::VirtualBox)

    xx = (rand() .- 0.5) .* 2 .*Vbinfo.σ
    yy = (rand() .- 0.5) .* (Vbinfo.Y_end - Vbinfo.Y_start) .+ (Vbinfo.Y_end + Vbinfo.Y_start) ./ 2
    zz =(rand() .- 0.5) .* (Vbinfo.Z_end - Vbinfo.Z_start) .+ (Vbinfo.Z_end + Vbinfo.Z_start) ./ 2

    return[xx,yy,zz]
end




function uᵢ(vec_points::Vector{Vector{Float64}}, ϵᵢ::Float64, xᵢ::Vector{Float64}, σ::Float64)
    map(x -> ϵᵢ .* fσ((x .- xᵢ)./σ ), vec_points)
end
"""
Compute the new position of all the Eddies. We consider only the convective velocity along x axis. If outside the Virtual Box, a new eddy is randomly generated inside the Virtual Box
"""
function convect_eddy(dt, Eddy, U₀, σ, Vbinfo)
    x_tmp = Eddy.xᵢ[1] + dt * U₀
    if x_tmp < σ
        Eddy.xᵢ = [x_tmp, Eddy.xᵢ[2], Eddy.xᵢ[3]]
    else
        Eddy.xᵢ = new_rand_position(Vbinfo)
        Eddy.ϵᵢ = rand((-1,1), 3)
    end
    return Eddy
end

"""
The velocity in the 3 directions is computed in each point provided in x
"""
function compute_uᵢₚ(x::Vector{Vector{Float64}}, dt::Float64, Eddies::Vector{SEM_EDDY}, U₀::Float64, Vbinfo::VirtualBox)
    contribution = zeros(length(x),3)
    for j = 1:1:length(Eddies)
        Eddies[j] = convect_eddy(dt, Eddies[j], U₀, Vbinfo.σ,Vbinfo)
        contribution[:,1] .+= uᵢ(x, Eddies[j].ϵᵢ[1], Eddies[j].xᵢ, Vbinfo.σ)
        contribution[:,2] .+= uᵢ(x, Eddies[j].ϵᵢ[2], Eddies[j].xᵢ, Vbinfo.σ)
        contribution[:,3] .+= uᵢ(x, Eddies[j].ϵᵢ[3], Eddies[j].xᵢ, Vbinfo.σ)

    end

    return sqrt(Vbinfo.V_b/(Vbinfo.σ^3)) ./ (Vbinfo.N)^0.5 .* contribution, Eddies
end



"""
Create a vector of points
"""
function create_vector_points(x, y, z)
    vector_points = Vector{Float64}[]
    for i = 1:1:length(x)
        for j = 1:1:length(y)
            for k = 1:1:length(z)
                push!(vector_points, [x[i], y[j], z[k]])
            end
        end
        
    end
    return vector_points
end


"""
Compute the acutual velocity and the turbulent kinetic energy. The convective velocity is just in the x direction
"""
function compute_U_k(q::Matrix{Float64}, A::Matrix{Float64}, U₀::Float64)
    U = A * q'
    U = U'
    U[1] = U[1] .+ U₀
    k = zeros(size(U)[1])
    for i = 1:1:size(U)[1]
        k[i] = 0.5.* (U[i,1].^2 .+ U[i,2].^2 .+ U[i,3].^2)

    end
    return U, k    
end