using LinearAlgebra

function fσ(x, shape_fun)
    shape_fun(x)
end

"""
Tent-like shape function. The domain is [-1,1]x[-1,1]x[-1,1].
It satisfy the normalization condition:
```math
int_{-1}^{1} fσ^2(x) dx = 1
```
"""
function tent_fun(x)
    if abs(x[1]) <= 1 && abs(x[2]) <= 1 && abs(x[3]) <= 1
        
        return sqrt(1.5) * (1 - abs(x[1])) * 
         sqrt(1.5) * (1 - abs(x[2])) *
         sqrt(1.5) * (1 - abs(x[3]))
    
    else
        return 0
    end
end



"""
Step function. The domain is [-1,1]x[-1,1]x[-1,1].
It satisfy the normalization condition:
```math
int_{-1}^{1} fσ^2(x) dx = 1
```
"""
function step_fun(x)
    if abs(x[1]) <= 1 && abs(x[2]) <= 1 && abs(x[3]) <= 1

        return 1/sqrt(8)
    else
        return 0
        

    end
end



"""
Truncated Gaussian function. The domain is [-1,1]x[-1,1]x[-1,1].
It satisfy the normalization condition:
```math
int_{-1}^{1} fσ^2(x) dx = 1
```
"""
function trunc_gauss_fun(x)
    if abs(x[1]) <= 1 && abs(x[2]) <= 1 && abs(x[3]) <= 1 #
        Cm = 1.301^3 #0.7511 #0.574872
        return Cm*exp(-9*(x[1])^2/ 2)*exp(-9*(x[2])^2/ 2)*exp(-9*(x[3])^2/ 2)
    else
        return 0        

    end
end





