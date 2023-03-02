"""
    cholesky_decomposition(R::Matrix{Float64})
    
Cholesky Decomposition of the Reynolds Stress Tensor.
"""
function cholesky_decomposition(R::Matrix{Float64})
    if length(R) == 4 #For 2D case
        error("The Reynolds stress tensor has to be a 3x3 matrix, 2x2 matrix not supported")

        # a11 = (R[1, 1])^0.5
        # a21 = (R[2, 1]) / a11
        # a22 = (R[2, 2] - a21)^0.5

        # A = [a11 0.0
        #     a21 a22]
    elseif length(R) == 9 #For 3D case
        a11 = (R[1, 1])^0.5
        a21 = (R[2, 1]) / a11
        a22 = (R[2, 2] - a21^2)^0.5
        a31 = (R[3, 1]) / a11
        a32 = (R[3, 2] - a21 * a31) / a22
        a33 = (R[3, 3] - a31^2 - a32^2)^0.5
        A = [a11 0.0 0.0
            a21 a22 0.0
            a31 a32 a33]
    else
        error("The Reynolds stress tensor has to be a 3x3 matrix")
    end

    return A
end
