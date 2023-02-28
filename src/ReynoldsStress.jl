using Interpolations
using DataFrames,XLSX
using SyntheticEddyMethod

"""
Function called by the user, where the Re_file_info is path of the .xlsx file with the 
data of the Reynolds Stress
File 
Z Y UU VV WW UV UW VW
"""
 function get_reynolds_stress_from_file(Re_file_info::String)
    Reinfo = DataFrame(XLSX.readtable(Re_file_info, "Sheet1")...)
    yy,zz = get_unique_coordinates_from_file(Reinfo)
    Re_tensor = Reynolds_stress_tensor(length(yy),length(zz))
    
    for i in eachindex(Reinfo.Z)
        idy = findall(x-> x == Reinfo.Y[i], yy)[1]
        idz = findall(x-> x == Reinfo.Z[i], zz)[1]
        Re_tensor.UU[idy,idz] = Reinfo.UU[i]
        Re_tensor.VV[idy,idz] = Reinfo.VV[i]
        Re_tensor.WW[idy,idz] = Reinfo.WW[i]
        Re_tensor.UV[idy,idz] = Reinfo.UV[i]
        Re_tensor.UW[idy,idz] = Reinfo.UW[i]
        Re_tensor.VW[idy,idz] = Reinfo.VW[i]
    end
    Reynolds_stress_interpolator(yy, zz, Re_tensor)
 end




function get_unique_coordinates_from_file(Reinfo::DataFrame)
    col_names = names(Reinfo)

    if isempty(findall(x-> x == "Y", col_names))
        error("In Reynolds Stress File, column Y not found")
    end

    if isempty(findall(x-> x == "Z", col_names))
        error("In Reynolds Stress File, column Z not found")
    end

    yy = convert(Array{Float64,1}, unique(Reinfo.Y)) 
    zz = convert(Array{Float64,1}, unique(Reinfo.Z)) 
    return yy, zz
end

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
end


function Reynolds_stress_tensor(a,b)
    Reynolds_stress_tensor(zeros(a,b),zeros(a,b),zeros(a,b),zeros(a,b),zeros(a,b),zeros(a,b))  
end

function Reynolds_stress_interpolator(yy::Vector{Float64},zz::Vector{Float64}, S::Reynolds_stress_tensor)
    Reynolds_stress_interpolator(LinearInterpolation((yy, zz), S.UU), 
    LinearInterpolation((yy, zz), S.VV), LinearInterpolation((yy, zz), S.WW), 
    LinearInterpolation((yy, zz), S.UV), LinearInterpolation((yy, zz), S.UW),
    LinearInterpolation((yy, zz), S.VW))
end


"""
Create vec_A = Vector{Matrix}, length(vec_A) == length(vec_points)
Each element is the Reynolds stress at the specific point
"""
function Reynolds_stress_points(vec_points::Vector{Vector{Float64}}, A::Reynolds_stress_interpolator)
    # vec_A = map(x ->Reynolds_stress_point(x, A), vec_points)
    npoints= length(vec_points)
    vec_A = Vector{Matrix}(undef,npoints)
    for i = 1:1:npoints
        vec_A[i] = Reynolds_stress_point(vec_points[i], A)
    end
    return vec_A
end

function Reynolds_stress_points(vec_points::Vector{Vector{Float64}}, A::Matrix{Float64})
    # vec_A = map(x -> A, vec_points )
    npoints= length(vec_points)
    vec_A = Vector{Matrix}(undef,npoints)
    for i = 1:1:npoints
        vec_A[i] = A
    end
    return vec_A
end

"""
Compute the interpolation of the point pp. It extracts the second and third coordinate.
"""
function Reynolds_stress_point(pp::Vector{Float64}, Re::Reynolds_stress_interpolator)
    y = pp[2]
    z = pp[3]
    Stress_Mat = [Re.UU(y,z) Re.UV(y,z) Re.UW(y,z); 
    Re.UV(y,z) Re.VV(y,z) Re.VW(y,z); 
    Re.UW(y,z) Re.VW(y,z) Re.WW(y,z)]
    return Stress_Mat
end
z = [[0.0497595595233471, -0.005560114405394889, -0.00035034087745518624], [0.049506494019851965, -0.005596091852882159, -0.0003933124821482676], [0.049253428516356834, -0.005632069300369428, -0.0004362840868413499], [0.0490003630128617, -0.0056680467478567, -0.0004792556915344221], [0.048747297509366794, -0.005704024195343966, -0.0005222272962275032], [1.0484942320058717, -0.005740001642831242, -0.0005651989009205832], [1.0482411665023768, -0.005775979090318509, -0.000608170505613663], [1.0479881009988816, -0.005811956537805783, -0.0006511421103067365], [1.0477350354953865, -0.005847933985293048, -0.0006941137149998228], [1.0474819699918916, -0.005883911432780317, -0.000737085319692901]]
mapreduce(permutedims, vcat, z)

