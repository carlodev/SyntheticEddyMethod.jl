using SyntheticEddyMethod
using Test

using Statistics
using DataFrames, XLSX, Plots


## Power Spectral Density Analysis
TI_vec = 0.1:0.1:0.5
k = 0.1:1000
E = (k).^(-5/3) .*100 #multiplied by 100 for shifting the curve in the top part

σ = 0.1 #eddy dimensions, the same in all the directions
b = 5.0
a = 0.0
#Defining the Virtual Box domain
x = collect(-σ:0.1:+σ) 
y = collect(a:0.1:b)
z = collect(a:0.1:b)


N_restart = 20
freqs = 0.0   
Nt = 2000
PSD = 0.0
vector_points = [[0.0, b/2, b/2]]




Vboxinfo = VirtualBox(y,z,σ;shape_fun=:trunc_gauss)


U₀ = 1.0 #Convective Velocity
dt = 0.001

#It can take up to 10 minutes
for TI in TI_vec
    A, Eddies = initialize_eddies(U₀, TI, Vboxinfo)

    for i=1:1:N_restart
        q = zeros(Nt, 3)
        for i = 1:1:Nt
            q[i,:] = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]
        
        end


        U, Ek =  compute_U_k(q, A, U₀)
        PSD_tmp, freqs = fft_from_signal(Ek, dt)
        PSD = PSD .+ PSD_tmp./N_restart

    end
        PSD_data = DataFrame([PSD, freqs], [:PSD, :freqs])
        XLSX.writetable("test/psd_results_$TI.xlsx", "$TI" => PSD_data)
end

#Random Signal
N_rand = 1000
PSD_rand_tot = 0.0
freqs_rand = 0.0
for i = 1:1:N_rand
    rand_signal = randn(3000)
    PSD_rand, freqs_rand = fft_from_signal(3/2 .* rand_signal.^2 ,dt)
    PSD_rand_tot = PSD_rand_tot .+ 1/N_rand .*PSD_rand
end




PSD_data = DataFrame[]

for i = eachindex(TI_vec)
    TI = TI_vec[i]
    filename = "test/psd_results_$TI.xlsx"
    df_tmp = DataFrame(XLSX.readtable(filename, "$TI"))
    push!(PSD_data, df_tmp)
end


Plots.plot(xaxis=:log, yaxis=:log, xlim = [0.5, 1e3], ylims =[1e-7, 1e2], xlabel="k", ylabel="E(k)", legend=:bottomleft, xticks=[1,10,100,1000])
for i = eachindex(TI_vec)
    TI = TI_vec[i]
    Plots.plot!(PSD_data[i].freqs, PSD_data[i].PSD, label = "SEM - TI = $TI")

end

plot!(freqs_rand, PSD_rand_tot, label = "RAND")
plot!(k, E, linestyle=:dash, label = "E(k)∝k^-5/3")


#Plots.plot!(freqs, PSD, label = "SEM - TI = $TI")
#Plots.plot!(freqs_mean, PSD_mean, label = "SEM mean")
Plots.savefig("SEM_vs_RAND_gaus.png")