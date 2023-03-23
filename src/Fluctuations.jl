"""
    compute_uSEM(vec_points::Vector{Vector{Float64}}, Eddies::Vector{SemEddy}, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator})

It computes the velocity fluctuations using the SEM. In order it computes
- Each j-eddy is convected of a distance dt⋅U₀ in the x direction
- The distance between each point and the centre of the eddy x - xⱼ
- It is normalised using σ for each direction
- qσ using the shape function and taking into account the intensity ϵᵢ in each direction
-  Reynolds stress and the Cholesky decomposition
- total contribution of the j-eddy

At the end the total contribution is rescaled by a factor B
"""
function compute_uSEM(point::Vector{Float64}, Eddies::Vector{SemEddy}, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator})
    vec_rk = map(E -> compute_rk(point, E.xᵢ), Eddies)
    rk_σ = map((rk, E) -> rk ./ E.σ, vec_rk, Eddies)
    qσ = map((x, E) -> fσ(x, Vbinfo.shape_fun) .* E.ϵᵢ, rk_σ, Eddies)
    Ap = Reynolds_stress_points(point, Re)
    uk = map(x -> Ap * x, qσ)
    contribution = sum(uk)

    σ_mean = Vbinfo.σ[1] * Vbinfo.σ[2] * Vbinfo.σ[3]

    B = sqrt(Vbinfo.V_b / (σ_mean)) ./ (Vbinfo.N)^0.5

    return B .* contribution, Eddies
end

function compute_uSEM(vec_points::Vector{Vector{Float64}}, Eddies::Vector{SemEddy}, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator})
    return map(x-> compute_uSEM(x, Eddies, Vbinfo, Re), vec_points)
end

"""
    compute_fluct(vec_points::Vector{Vector{Float64}}, dt::Float64, Eddies::Vector{SemEddy}, U₀::Float64, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator}; DFSEM = false)

Compute the velocity fluctuations accordingly to the Reynolds Stress `Re`. It can be selected the DFSEM that impose also the divergence free condition.
"""
function compute_fluct(point::Vector{Float64}, T::Float64, Eddies::Vector{SemEddy}, U₀::Float64, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator}; DFSEM=false)

    #Convect Eddies
    if T > Vbinfo.T
        @info "convecting eddies"
        dt = T-Vbinfo.T
        Eddies = map(ej -> convect_eddy(dt, ej, U₀, Vbinfo), Eddies)
        Vbinfo.T = T
    end

    if DFSEM
        U, Eddies = compute_uDFSEM(point, Eddies, Vbinfo, Re)
    else
        U, Eddies = compute_uSEM(point, Eddies, Vbinfo, Re)
    end

    # Add U₀, convective velocity, to the component in the x direction. Save the results back to U
    u_ =[U[1] + U₀, U[2], U[3]] 

    return u_
end

function compute_fluct(vec_points::Vector{Vector{Float64}}, dt::Float64, Eddies::Vector{SemEddy}, U₀::Float64, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator}; DFSEM=false)
    return map(x-> compute_fluct(x, dt, Eddies, U₀, Vbinfo, Re; DFSEM=DFSEM), vec_points)
end

    ### DFSEM




   #Compute a vector for the left term of the DFSEM
function compute_RL(Re::Matrix{Float64})
    eig_vals, eig_vec = eigen(Re)
    ReL = eig_vec.*eig_vals'
    #compute k 
    k = compute_kp(eig_vals)
    map!(ev -> verify_real_sqrt(k, ev), eig_vals, eig_vals)

    #Obtaining a Vector
    ReL =  collect((sqrt.(2 .* (k .- eig_vals))' * eig_vec')')
    return ReL
end 
#Compute a vector for the left term of the DFSEM
# function compute_RL(Re::Matrix{Float64})
#     eig_vals, ReLG = eigen(Re)
#     Re
#     k = compute_kp(Re)

#     map!(ev -> verify_real_sqrt(k, ev), eig_vals, eig_vals)

#     #Obtaining a Vector
#     ReL = ReLG * sqrt.(2 .* (k .- eig_vals))
#     return ReL
# end

function verify_real_sqrt(k::Real, ev::Real)
    if k - ev < 0
        @warn "Cannot reproduce the exact Reynolds stress for this point: k - ev <0, reducing ev from $ev tp $(0.5.*k)"
        ev = 0.5 .* k
    end

    return ev
end

#Compute the left side vector of the cross product in DFSEM
function compute_α(RL::Vector{Float64}, epsi::Vector)
    alpha = RL .* epsi
    return alpha
end

#Compute turbulent kinetic energy from the Reynols stress extracting the trace
function compute_kp(Re::Matrix{Float64})
    k = 0.5 * tr(Re)
    return k
end

#Compute turbulent kinetic energy from the Reynols stress extracting the trace
function compute_kp(eig_val::Vector{Float64})
    k = 0.5 * sum(eig_val)
    return k
end


#Compute the distance between a point and the centre of the eddy
function compute_rk(x::Vector{Float64}, xk::Vector{Float64})
    vec_rk = x - xk

    return vec_rk
end


"""
    compute_uDFSEM(vec_points::Vector{Vector{Float64}}, Eddies::Vector{SemEddy}, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator})

    It computes the velocity fluctuations using the DFSEM. In order it computes
- Each j-eddy is convected of a distance dt⋅U₀ in the x direction
- The distance between each point and the centre of the eddy x - xⱼ
- It is normalised using σ for each direction, takes the norm
- qσ using the specifically designed shape function
- eigenvalues of Reynolds stress in principal axes and its eigenvectors
- a rotation from Local (principal axes) to Global coordinates system, based on eigenvectors of Reynolds stress tensor and the eddy intensity ϵᵢ
- total contribution of the j-eddy

At the end the total contribution is rescaled by a factor B
"""
function compute_uDFSEM(point::Vector{Float64}, Eddies::Vector{SemEddy}, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator})
    # @assert Vbinfo.σ[1] ==  Vbinfo.σ[2] 
    # @assert  Vbinfo.σ[1] == Vbinfo.σ[3]
    if Vbinfo.shape_fun != DFSEM_fun
        Vbinfo.shape_fun = DFSEM_fun
        println("Set the shape functions a DFSEM_fun")
    end
    # @assert Vbinfo.shape_fun == DFSEM_fun
    @assert typeof(Re) == Matrix{Float64} #Not supported point defined reynolds for DFSEM

    # σ = Vbinfo.σ[1]
    #Using improperly the create_vector_points function, just to intialize the contribution of the eddies
    vec_rk = map(E -> compute_rk(point, E.xᵢ), Eddies)
    #compute rk/σ
    rk_σ = map(x -> norm(x ./ Vbinfo.σ), vec_rk)

    #Compute qσ
    qσ = map(x -> DFSEM_fun(x), rk_σ)

    r_vec = map((qσk, rk_σk, vec_rkk) -> qσk / (rk_σk^3) .* (vec_rkk ./ Vbinfo.σ), qσ, rk_σ, vec_rk)


    ReL = compute_RL(Re)

    l_vec = map(E -> compute_α(ReL, E.ϵᵢ), Eddies)
    uk = map((r_x, l_x) -> cross(r_x, l_x), r_vec, l_vec)
    contribution = sum(uk)

    σ_mean3 = Vbinfo.σ[1] * Vbinfo.σ[2] * Vbinfo.σ[3]
    B = (15 * Vbinfo.V_b / (16 * pi * σ_mean3))^0.5

    return B .* (1 ./ (Vbinfo.N)^0.5) .* contribution, Eddies
end

function compute_uDFSEM(vec_points::Vector{Vector{Float64}}, dt::Float64, Eddies::Vector{SemEddy}, U₀::Float64, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator}; DFSEM=false)
    return map(x-> compute_uDFSEM(x, dt, Eddies, U₀, Vbinfo, Re; DFSEM=DFSEM), vec_points)
end