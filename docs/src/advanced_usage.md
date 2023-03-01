# Advanced usage
In this page are provided some example on how the user can customize the eddy generation.

Defining the Virtual Box in `y` and `z` directions.
```julia
using SyntheticEddyMethod

σ = 0.1
b = 2.0
a = 0.0

y = collect(a:0.1:b)
z = collect(a:0.1:b)
```

It is possible also to specify the `x` direction dimensions, a different shape funcion and different `σ` in the three different directions.
The shape coded can be seen in the `Shapefunctions.jl` file.

```julia
x = collect(-1:0.1:1)
σx = 0.1
σy = 0.05
σz = 0.07
σ = [σx,σy,σz]
Vboxinfo = VirtualBox(x,y,z,σ; shape_fun = step_fun)
N = Vboxinfo.N
Eddies = initialize_eddies(N, σ, Vboxinfo)
```

### Define custom Reynolds Stress Tensor

The user can also use a custom Reynolds stress just by writing a `3x3` matrix.
```julia
A = [1.87e-5 -8.6e-8  -8.2e-7; -8.6e-8 5.27e-8 6.8e-8; -8.2e-7 4.9e-9 2.64e-6]
```

### Import the Reynolds Stress Tensor from file
Or use a database where the the Reynolds Stress is defined pointwise. 
```julia
reynolds_stress_file = joinpath(@__DIR__,"..","src","Data","Re_ch.xlsx")
A_from_file = get_reynolds_stress_from_file(reynolds_stress_file)
```

An curious user notice that in this last case the Reynolds Stress is not a matrix, but it is an interpolator object. It depends on the point location where the user want to compute the fluctuations.

```julia
vector_points = [[0.0, 1.0, 2.5]]
dt = 0.01
U₀ = 1.0
compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo,A )
```


# Analyze the signal
A simple case is used.
```julia
using SyntheticEddyMethod
using Statistics

σ = 0.1 #eddy dimensions, the same in all the directions
b = 5.0
a = 0.0


#Defining the Virtual Box domain
x = -σ:0.1:+σ 
y = collect(a:0.1:b)
z = collect(a:0.1:b)


Vboxinfo = VirtualBox(y,z,σ)
N = Vboxinfo.N
dt = 0.01

U₀ = 1.0 #Convective Velocity
TI = 0.01 #turbulence intensity

A, Eddies = initialize_eddies(U₀, TI, Vboxinfo)


#Computing the velocity in the middle of the VirtualBox domain
vector_points = [[0.0, b/2, b/2]]

#Defining how many time interval
Nt = 1000
U = zeros(Nt, 3)


for i = 1:1:Nt
    U[i,:] = compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo, A)[1]
end

#The deviation standard should approach the turbulence intensity
Statistics.std(U[:,1])
Statistics.std(U[:,2])
Statistics.std(U[:,3])

```
## Spectral Analysis
```julia
using SyntheticEddyMethod
using DataFrames, XLSX, Plots
#Turbulence intensity to test
TI_vec = 0.01:0.01:0.05

k = 0.1:1000
E = (k).^(-5/3)*0.01 #multiplied by 100 for shifting the curve in the top part

N_restart = 20
Nt = 2000
dt = 0.001

#It can take up to 10 minutes
for TI in TI_vec
    A, Eddies = initialize_eddies(U₀, TI, Vboxinfo)
    PSD = 0.0
    freqs = 0.0   

    for i=1:1:N_restart
         q = Vector{Float64}[]
        for j = 1:1:Nt
            qi = compute_fluct(vector_points, dt, Eddies, U₀, Vboxinfo, A)[1]
            push!(q,qi)
        end
        println(i)

        Ek =  compute_Ek(q, U₀)
        PSD_tmp, freqs = fft_from_signal(Ek, dt)
        PSD = PSD .+ PSD_tmp ./N_restart

    end

PSD_data = DataFrame([PSD, freqs], [:PSD, :freqs])
        XLSX.writetable("test/psd_results_$TI.xlsx", "$TI" => PSD_data)
end

#Random Signal
N_rand = 1000
PSD_rand_tot = 0.0
freqs_rand = 0.0
for i = 1:1:N_rand
    rand_signal = randn(3000).*(TI)
    PSD_rand, freqs_rand = fft_from_signal(3/2 .* rand_signal.^2 ,dt)
    PSD_rand_tot = PSD_rand_tot .+ 1/N_rand .*PSD_rand
end




PSD_data = DataFrame[]

for i = eachindex(TI_vec)
    TI = TI_vec[i]
    filename = "test/psd_results_$TI.xlsx"
    df_tmp = DataFrame(XLSX.readtable(filename, "$TI")...)
    push!(PSD_data, df_tmp)
end


Plots.plot(xaxis=:log, yaxis=:log, xlim = [0.5, 1e3], ylims =[1e-10, 1], xlabel="k", ylabel="E(k)", legend=:bottomleft, xticks=[1,10,100,1000])
for i = eachindex(TI_vec)
    TI = TI_vec[i]
    Plots.plot!(PSD_data[i].freqs, PSD_data[i].PSD, label = "SEM - TI = $TI")

end

plot!(freqs_rand, PSD_rand_tot, label = "RAND")
plot!(k, E, linestyle=:dash, label = "E(k)∝k^-5/3")


```


