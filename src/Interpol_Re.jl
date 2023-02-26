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
        idy = Int(findall(x-> x == Reinfo.Y[i], yy)[1])
        idz = Int(findall(x-> x == Reinfo.Z[i], zz)[1])
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



refile ="src/Data/Re_ch.xlsx"
A = get_reynolds_stress_from_file(refile)





function Reynolds_stress_point(pp::Vector{Float64}, Re::Reynolds_stress_interpolator)
    y = pp[2]
    z = pp[3]
    Stress_Mat = [Re.UU(y,z) Re.UV(y,z) Re.UW(y,z); 
    Re.UV(y,z) Re.VV(y,z) Re.VW(y,z); 
    Re.UW(y,z) Re.VW(y,z) Re.WW(y,z)]
    return Stress_Mat
end

x = collect(1:1:5)
y = collect(0.01:0.01:1.5)
z = collect(0.1:0.1:6)
vv = create_vector_points(x,y,z)

map(x ->Reynolds_stress_point(x, A), vv )
