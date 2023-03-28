#Spectral Tools

"""
    fft_from_signal(q::Vector{Float64},dt::Float64)

It performs the Fast Fourier transformation of the signal `q`.

# Examples
```julia-repl
# We create a random signal
julia> q = rand(1000)
1000-element Vector{Float64}:
 0.7478100061738703
 0.16303784021483914
 0.3628099523805166
 0.2608705413573018
 0.22498693840252693
 0.9428769056802648
 0.3579787963511998
 0.7224804012744513
 ⋮
 0.9946122085650284
 0.9183278219776294
 0.19908223170935013
 0.7708756376182156
 0.5007439030983675
 0.8464298265686861
 0.8700496116560734
 0.3670620181980453
julia> psd, freq = fft_from_signal(q,0.1)
([0.11064559574015868, 0.005389219813554343, 0.21634910063286397, 0.12896616563010807, 0.1578838809521174, 0.051975560320304294, 0.10663607461615124, 0.1119872950890813, 0.009687265867758152, 0.1469444321528952  …  0.08051345232908684, 0.1481450312902375, 0.018682727272637666, 0.19356783851443643, 0.006863014095143317, 0.009106381063590546, 0.01426246439106842, 0.25854085870704596, 0.08248723638591579, 0.18281921062421522], [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1  …  4.9, 4.91, 4.92, 4.93, 4.94, 4.95, 4.96, 4.97, 4.98, 4.99])
julia> length(psd)
499
julia> length(freq)
499
```

"""
function fft_from_signal(q::Vector{Float64},dt)

    nt=length(q)
    fhat=fft(q)
    
    PSD = fhat.*conj(fhat)/(nt)
    PSD = real(fftshift(PSD))
    freqs = fftshift(fftfreq(nt,1/dt))
    idx = findall(x -> x>0, freqs)
    
    return PSD[idx], freqs[idx]
end



# VecPoints = Union{Float64, Vector{Float64}}

"""
    create_vector_points(x, y, z)

Create a vector of points. Useful for testing purposes.

# Examples
```julia-repl
julia> create_vector_points([1.0], [2.0, 3.0], [1.5, 3.5, 4.2])
6-element Vector{Vector{Float64}}:
 [1.0, 2.0, 1.5]
 [1.0, 2.0, 3.5]
 [1.0, 2.0, 4.2]
 [1.0, 3.0, 1.5]
 [1.0, 3.0, 3.5]
 [1.0, 3.0, 4.2]
```
"""
function create_vector_points(x, y, z)
    vector_points = Vector{Float64}[]
    for i = 1:1:length(x)
        for j = 1:1:length(y)
            for k = 1:1:length(z)
                push!(vector_points, [x[i], y[j], z[k]])
            end
        end
        
    end
    return vector_points
end


"""
    compute_Ek(U::Vector{Vector{Float64}}, U₀::Float64)

Compute the turbulent kinetic energy. The convective speed ` U₀` is subtracted from the `x` component of the speed.
# Examples
```julia-repl
julia> compute_Ek([1.2, 0.1, 0.3], 1.0)
0.06999999999999999
```
"""
function compute_Ek(U::Vector{Float64}, U₀::Float64)
    U = [U[1] - U₀, U[2], U[3]]
    Ek = 0.5 .* sum(U.^2)
    return Ek
end

#compute_Ek for Vector{Vector{Float64}} of points
function compute_Ek(U::Vector{Vector{Float64}}, U₀::Float64)
    return map(x->compute_Ek(x, U₀), U )
end