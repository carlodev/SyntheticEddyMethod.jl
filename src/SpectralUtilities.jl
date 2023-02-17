#Spectral Tools

function fft_from_signal(q,dt)
    nt=length(q)
    fhat=fft(q)
    
    PSD = fhat.*conj(fhat)/(nt)
    PSD = real(fftshift(PSD))
    freqs = fftshift(fftfreq(nt,1/dt))
    idx = findall(x -> x>0, freqs)

    return PSD[idx], freqs[idx]
end