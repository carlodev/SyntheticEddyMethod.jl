using LinearAlgebra

function fσ(x, shape_fun::Symbol)
    if shape_fun == :tent
        return tent_fun(x)
    elseif shape_fun == :step
        return step_fun(x)
    elseif shape_fun == :trunc_gauss
        return trunc_gauss_fun(x)
    end
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



# xtr = -1.3:0.001:1.3
# ytr  = -1.3:0.001:1.3
# ztr = ytr
# pp = Any[]
# for p in zip(xtr,ytr,ztr)
#     pv = [p[1], p[2], p[3]]
#     push!(pp,pv)
# end

# using Plots

# plot(xtr, (tent_fun.(pp)), label = "tent function", linewidth=2)
# savefig("tent_fun.pdf")
# plot(xtr, (step_fun.(pp)), legend=:bottom, label = "step function", linewidth=2)
# savefig("step_fun.pdf")
# plot(xtr, (trunc_gauss_fun.(pp)), label = "gaussian function", linewidth=2)
# savefig("trun_gauss_fun.pdf")


