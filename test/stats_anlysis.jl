using MCMCChains
using StatsPlots
using FileIO, JLD2
using DataFrames
using Statistics


Stat = FileIO.load("StatF2.jld2","Stat")
U = FileIO.load("UF2.jld2","U")




# construct a Chains object
chn = Chains(Stat)
chn_stat =Chains([Stat[:,1];Stat[:,3];Stat[:,2]], ["fluctuations"])
chn_U = Chains(U,["u'", "v'", "w'"])
Statistics.std([Stat[:,1];Stat[:,3];Stat[:,2]])

# Visualize the mean value of each components
StatsPlots.plot(chn, seriestype = :mixeddensity; size=(840, 600))

Plots.plot(chn_stat, seriestype = :mixeddensity, linewidth=2; size=(840, 600))
savefig("SEM_Fluctuations.pdf")

# Visualize the mean convergence for each componenents
meanplot(chn_U, linewidth = 2, size=(840, 600))
savefig("SEM_mean.pdf")



# Autocorrelation
using StatsBase
x_val_corr = collect(1:1:44) .*0.01

autocorrU = []
labU = ["u'", "v'", "w'"]
markerU = [:solid,:dashdot, :dash]
for i = 1:3
    push!(autocorrU,autocor(U[:,i]; demean=true))
end
autocorrU


Plots.plot(xlabel = "τ [s]", ylabel = "autocorrelation")
for i = 1:3
    autocorr_plt = Plots.plot!(x_val_corr, autocorrU[i], label = labU[i], linewidth=2, linestyle = markerU[i])
end
Plots.plot!(legend=:right)
# savefig("SEM_Autocorrelation.svg")

#Computing integral time
using Trapz

autocorrU[1]
trapz(x_val_corr, autocorrU[3])