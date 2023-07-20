"""
    get_reynolds_stress_from_file(Re_file_info::String)

Function called by the user, where the Re_file_info is path of the .xlsx file with the data of the Reynolds Stress.\\
File column example:

| Z | Y | UU | VV | WW | UV | UW | VW |
"""
 function get_reynolds_stress_from_file(Re_file_info::String; dims = (:Z,:Y))
    dims =  convert_symb2tupl(dims)
    Reinfo = DataFrame(XLSX.readtable(Re_file_info, "Sheet1"))
    yy,zz = get_unique_coordinates_from_file(Reinfo, dims)
    
    Re_tensor = Reynolds_stress_tensor(length(yy),length(zz))
    expr = :($Reinfo.$(dims[1]))
    Reidx = eval(expr)

    for i in eachindex(Reidx)

        if length(yy)>1
            idy = findall(x-> x == Reinfo.Y[i], yy)[1]
        else
            idy = 1

        end
        
        if length(zz)>1
            idz = findall(x-> x == Reinfo.Z[i], zz)[1]
        else
            idz = 1
        end

        Re_tensor.UU[idy,idz] = Reinfo.UU[i]
        Re_tensor.VV[idy,idz] = Reinfo.VV[i]
        Re_tensor.WW[idy,idz] = Reinfo.WW[i]
        Re_tensor.UV[idy,idz] = Reinfo.UV[i]
        Re_tensor.UW[idy,idz] = Reinfo.UW[i]
        Re_tensor.VW[idy,idz] = Reinfo.VW[i]
    end
    
    Reynolds_stress_interpolator(yy, zz, Re_tensor, dims)
end

function convert_symb2tupl(a::Symbol)
    a = (a, )
    return a
end

function convert_symb2tupl(a::Tuple)
    return a
end


function get_unique_coordinates_from_file(Reinfo::DataFrame, dims::Tuple)
    col_names = names(Reinfo)
    valid_cols = (:Y, :Z)
    for vc in valid_cols
        if isempty(findall(x-> x == "$vc", col_names)) && !isempty(findall(x-> x == vc, dims))
            error("In Reynolds Stress File, column $vc not found")
        end
    end

    if length(dims) == 1
        idx = dims[1]
        if idx == :Y
            yy = convert(Array{Float64,1}, unique(Reinfo.Y)) 
            zz = [0.0]
        elseif idx == :Z
            yy = [0.0]
            zz = convert(Array{Float64,1}, unique(Reinfo.Z)) 
        end
    elseif length(dims) == 2
        yy = convert(Array{Float64,1}, unique(Reinfo.Y)) 
        zz = convert(Array{Float64,1}, unique(Reinfo.Z)) 
    end
    
    return yy, zz
end


"""
    Reynolds_stress_tensor

Struct hosting the Reynolds stress tensor informations
"""
mutable struct Reynolds_stress_tensor
    UU::Matrix{Float64}
    VV::Matrix{Float64}
    WW::Matrix{Float64}
    UV::Matrix{Float64}
    UW::Matrix{Float64}
    VW::Matrix{Float64}
end

mutable struct Reynolds_stress_interpolator
    UU
    VV
    WW
    UV
    UW
    VW
    dims::Tuple
end


function Reynolds_stress_tensor(a::Int,b::Int)
    Reynolds_stress_tensor(zeros(a,b),zeros(a,b),zeros(a,b),zeros(a,b),zeros(a,b),zeros(a,b))  
end

function Reynolds_stress_interpolator(yy::Vector{Float64},zz::Vector{Float64}, S::Reynolds_stress_tensor, dims::Tuple)
    if length(dims) == 2
        inter_points = (yy,zz)
        return Reynolds_stress_interpolator(LinearInterpolation(inter_points, S.UU), 
        LinearInterpolation(inter_points, S.VV), LinearInterpolation(inter_points, S.WW), 
        LinearInterpolation(inter_points, S.UV), LinearInterpolation(inter_points, S.UW),
        LinearInterpolation(inter_points, S.VW), dims)

    elseif length(yy) > 1
        inter_points = yy
        return Reynolds_stress_interpolator(LinearInterpolation(inter_points, S.UU[:,1]), 
        LinearInterpolation(inter_points, S.VV[:,1]), LinearInterpolation(inter_points, S.WW[:,1]), 
        LinearInterpolation(inter_points, S.UV[:,1]), LinearInterpolation(inter_points, S.UW[:,1]),
        LinearInterpolation(inter_points, S.VW[:,1]), dims)

    else
        inter_points = zz
        return Reynolds_stress_interpolator(LinearInterpolation(inter_points, S.UU[1,:]), 
        LinearInterpolation(inter_points, S.VV[1,:]), LinearInterpolation(inter_points, S.WW[1,:]), 
        LinearInterpolation(inter_points, S.UV[1,:]), LinearInterpolation(inter_points, S.UW[1,:]),
        LinearInterpolation(inter_points, S.VW[1,:]), dims)
    end

end


"""
    Reynolds_stress_points(vec_points::Vector{Vector{Float64}}, A::Reynolds_stress_interpolator)

Create vec_A of type Vector{Matrix}, where length(vec_A) == length(vec_points)
Each element is the cholesky decomposition of the Reynolds stress at the specific point
"""
function Reynolds_stress_points(point::Vector{Float64}, Re_interp::Reynolds_stress_interpolator)
    # vec_A = map(x ->Reynolds_stress_point(x, A), vec_points)
    # npoints= length(vec_points)
    # vec_A = Vector{Matrix}(undef,npoints)
    # for i = 1:1:npoints
    #     Re_loc = Reynolds_stress_point(vec_points[i], Re_interp)
    #     vec_A[i] = cholesky_decomposition(Re_loc)
    # end

    Re_loc = Reynolds_stress_point(point, Re_interp)
    A = cholesky_decomposition(Re_loc)
    return A
end

"""
    Reynolds_stress_points(vec_points::Vector{Vector{Float64}}, A::Matrix{Float64})
    
Create vec_A of type Vector{Matrix}, where length(vec_A) == length(vec_points)
Each element is the cholesky decomposition of the Reynolds stress at the specific point.
"""
function Reynolds_stress_points(point::Vector{Float64}, Re::Matrix{Float64})

    # # vec_A = map(x -> A, vec_points )
    # npoints= length(vec_points)
    # vec_A = Vector{Matrix}(undef,npoints)
    # for i = 1:1:npoints
    #     vec_A[i] = cholesky_decomposition(Re)
    # end
    A = cholesky_decomposition(Re)
    return A
end

"""
    Reynolds_stress_point(point::Vector{Float64}, Re::Reynolds_stress_interpolator)

It gives the Reynolds stress for a specific point. It performs an interpolation using the information in Re.
"""
function Reynolds_stress_point(pp::Vector{Float64}, Re::Reynolds_stress_interpolator)
    y = pp[2]
    z = pp[3]
    if length(Re.dims)==1
        if Re.dims[1] == :Y
            eval_point = y
        elseif Re.dims[1] == :Z
            eval_point = z
        end
        Stress_Mat = [Re.UU(eval_point) Re.UV(eval_point) Re.UW(eval_point); 
        Re.UV(eval_point) Re.VV(eval_point) Re.VW(eval_point); 
        Re.UW(eval_point) Re.VW(eval_point) Re.WW(eval_point)]
    else
        Stress_Mat = [Re.UU(y,z) Re.UV(y,z) Re.UW(y,z); 
        Re.UV(y,z) Re.VV(y,z) Re.VW(y,z); 
        Re.UW(y,z) Re.VW(y,z) Re.WW(y,z)]
    end
  
    return Stress_Mat
end
