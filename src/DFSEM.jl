using LinearAlgebra
function DFSEM_fun(rs, B)
    if norm(rs) <= 1 

        return  B*(sin(pi*rs)^2)*rs
        # return 1 -rs^2
    else
        return 0
    end
end

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

function compute_α(RL::Vector{Float64}, epsi::Vector)
    alpha = RL .* epsi
    return alpha
end


function compute_kp(Re::Matrix{Float64})
    k = 0.5*tr(Re)
    return k
end

function compute_rk(x::Vector{Float64}, xk::Vector{Float64})
    vec_rk = x-xk
    
    return vec_rk
end



function compute_uDFSEM(vec_points::Vector{Vector{Float64}}, dt::Float64, Eddies::Vector{SemEddy}, U₀::Float64, Vbinfo::VirtualBox, Re::Union{Matrix,Reynolds_stress_interpolator})
    σ_mean3 =  Vbinfo.σ[1] *  Vbinfo.σ[2] * Vbinfo.σ[3] 
    B = (16*Vbinfo.V_b/(15*pi*σ_mean3))^0.5
    
    # B = sqrt(Vbinfo.V_b/(σ_mean3))

    @assert Vbinfo.σ[1] ==  Vbinfo.σ[2]
    @assert  Vbinfo.σ[1] == Vbinfo.σ[3] 
    @assert Vbinfo.shape_fun == DFSEM_fun
    @assert typeof(Re) == Matrix{Float64} #Not supported point defined reynolds for DFSEM
    σ = Vbinfo.σ[1]
    #Using improperly the create_vector_points function, just to intialize the contribution of the eddies
    contribution =  create_vector_points(zeros(length(vec_points)),0,0)

    for j = 1:1:length(Eddies)
        Eddies[j] = convect_eddy(dt, Eddies[j], U₀, Vbinfo)
        #get rk vector
        vec_rk = map(x -> compute_rk(x, Eddies[j].xᵢ), vec_points)
        #compute rk/σ
        rk_σ = norm.(vec_rk) ./ σ
        #Compute qσ
        qσ= map(x-> DFSEM_fun(x, B), rk_σ)
        
        r_vec = map((qσk,rk_σk,vec_rkk)-> qσk/(rk_σk^3).*vec_rkk./σ, qσ,rk_σ,vec_rk)
        
      
        ReL = compute_RL(Re)
        
        l_vec = compute_α(ReL, Eddies[j].ϵᵢ)
        uk =  map(r_x -> cross(r_x, l_vec), r_vec)
        contribution .+= uk

    end

    return (1 ./ (Vbinfo.N)^0.5) .* contribution, Eddies
end
