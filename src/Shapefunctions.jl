"""
Tent-like shape function. The domain is [-1,1]x[-1,1]x[-1,1].
It satisfy the normalization condition:
```math
int_{-1}^{1} fσ^2(x) dx = 1
```
"""
function fσ(x)
   
    if abs(x[1]) <= 1 && abs(x[2]) <= 1 && abs(x[3]) <= 1

        return sqrt(1.5) * (1 - abs(x[1])) *
        sqrt(1.5) * (1 - abs(x[2])) *
        sqrt(1.5) * (1 - abs(x[3]))
    else
        return 0
        

    end
end