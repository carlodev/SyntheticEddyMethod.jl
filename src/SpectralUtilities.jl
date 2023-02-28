#Spectral Tools

function fft_from_signal(q::Vector{Float64},dt)

    nt=length(q)
    fhat=fft(q)
    
    PSD = fhat.*conj(fhat)/(nt)
    PSD = real(fftshift(PSD))
    freqs = fftshift(fftfreq(nt,1/dt))
    idx = findall(x -> x>0, freqs)
    return PSD[idx], freqs[idx]
end

psd, fre = fft_from_signal(rand(5), 0.001)

psd
fre