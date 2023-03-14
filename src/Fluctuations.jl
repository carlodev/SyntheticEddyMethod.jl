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
function compute_uSEM(vec_points::Vector{Vector{Float64}}, Eddies::Vector{SemEddy}, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator})
    contribution =  create_vector_points(zeros(length(vec_points)),0,0)
    for j = 1:1:length(Eddies)
        vec_rk = map(x -> compute_rk(x, Eddies[j].xᵢ), vec_points)
        rk_σ = map(x-> x ./ Eddies[j].σ, vec_rk)
        qσ= map(x-> fσ(x, Vbinfo.shape_fun) .* Eddies[j].ϵᵢ, rk_σ)

        Ap = Reynolds_stress_points(vec_points, Re)
        uk = map((x,y) -> x * y, Ap, qσ)
        contribution .+= uk

    end

    σ_mean =  Vbinfo.σ[1] *  Vbinfo.σ[2] * Vbinfo.σ[3] 
    B = sqrt(Vbinfo.V_b/(σ_mean)) ./ (Vbinfo.N)^0.5
    return  B.* contribution, Eddies
end


"""
    compute_fluct(vec_points::Vector{Vector{Float64}}, dt::Float64, Eddies::Vector{SemEddy}, U₀::Float64, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator}; DFSEM = false)

Compute the velocity fluctuations accordingly to the Reynolds Stress `Re`. It can be selected the DFSEM that impose also the divergence free condition.
"""
function compute_fluct(vec_points::Vector{Vector{Float64}}, dt::Float64, Eddies::Vector{SemEddy}, U₀::Float64, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator}; DFSEM = false)
    
    #Convect Eddies
    Eddies = map(ej -> convect_eddy(dt, ej, U₀, Vbinfo), Eddies)

    if DFSEM
        U, Eddies = compute_uDFSEM(vec_points, Eddies, Vbinfo, Re)
    else
        U, Eddies = compute_uSEM(vec_points, Eddies, Vbinfo, Re)
        # u_fluct = compute_uᵢₚ(vec_points, dt, Eddies, U₀, Vbinfo)[1]
        # u_fluct_vec = [u_fluct[i,:] for i in axes(u_fluct,1)]
    
        # Ap = Reynolds_stress_points(vec_points, Re)
        # U = map((x,y) -> x * y, Ap, u_fluct_vec)
    
    end
    
    # Add U₀, convective velocity, to the component in the x direction. Save the results back to U
    u_ = map(x -> [x[1] + U₀, x[2], x[3]], U)
    
    return u_ 
end


### DFSEM

#Compute a vector for the left term of the DFSEM
function compute_RL(Re::Matrix{Float64})
    # eig_vecs = eigvecs(Re)
    # eig_vals = eigvals(Re)
    # M = I(3).*eig_vals
    # ReLG = Re * inv(M)
    eig_vals, ReLG = eigen(Re)
    k = compute_kp(Re)
    for ev in eig_vals @assert k - ev >0 end

    #Obtaining a Vector
    ReL = ReLG * sqrt.(2 .*(k .- eig_vals))
    return ReL
end

#Compute the left side vector of the cross product in DFSEM
function compute_α(RL::Vector{Float64}, epsi::Vector)
    alpha = RL .* epsi
    return alpha
end

#Compute turbulent kinetic energy from the Reynols stress extracting the trace
function compute_kp(Re::Matrix{Float64})
    k = 0.5*tr(Re)
    return k
end

#Compute the distance between a point and the centre of the eddy
function compute_rk(x::Vector{Float64}, xk::Vector{Float64})
    vec_rk = x-xk
    
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
function compute_uDFSEM(vec_points::Vector{Vector{Float64}}, Eddies::Vector{SemEddy}, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator})
     # @assert Vbinfo.σ[1] ==  Vbinfo.σ[2] 
    # @assert  Vbinfo.σ[1] == Vbinfo.σ[3]
    if Vbinfo.shape_fun != DFSEM_fun
        Vbinfo.shape_fun = DFSEM_fun
        println("Set the shape functions a DFSEM_fun")
    end
    # @assert Vbinfo.shape_fun == DFSEM_fun
    @assert typeof(Re) == Matrix{Float64} #Not supported point defined reynolds for DFSEM
    σ = Vbinfo.σ[1]
    #Using improperly the create_vector_points function, just to intialize the contribution of the eddies
    contribution =  create_vector_points(zeros(length(vec_points)),0,0)

    for j = 1:1:length(Eddies)
        #get rk vector
        vec_rk = map(x -> compute_rk(x, Eddies[j].xᵢ), vec_points)
        #compute rk/σ
        rk_σ = map(x-> norm(x ./ Vbinfo.σ), vec_rk)

        #Compute qσ
        qσ= map(x-> DFSEM_fun(x), rk_σ)
        
        r_vec = map((qσk,rk_σk,vec_rkk)-> qσk/(rk_σk^3).*(vec_rkk./Vbinfo.σ), qσ,rk_σ,vec_rk)
        
      
        ReL = compute_RL(Re)
        
        l_vec = compute_α(ReL, Eddies[j].ϵᵢ)
        uk =  map(r_x -> cross(r_x, l_vec), r_vec)
        contribution .+= uk

    end

    σ_mean3 =  Vbinfo.σ[1] *  Vbinfo.σ[2] * Vbinfo.σ[3] 
    B = (15*Vbinfo.V_b/(16*pi*σ_mean3))^0.5
    
    return B.*(1 ./ (Vbinfo.N)^0.5) .* contribution, Eddies
end
