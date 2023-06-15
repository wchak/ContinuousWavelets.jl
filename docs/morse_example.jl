using ContinuousWavelets, Plots, Wavelets, FFTW, Logging, LinearAlgebra

#=

The Generalized Morse Wavelets (GMWs) is a superfamily of analytic wavelets.

The GMWs in the frequency domain can be defined by

$$\Psi_{\beta, \gamma}(\omega) := \int \psi_{\beta, \gamma}(t) e^{-i \omega t} d{t} = H(\omega) \alpha_{\beta, \gamma} \omega^\beta e^{-\omega^\gamma} $$

where $\beta>0, \gamma \geq 1$ are two main parameters controlling the form of wavelet, $$\alpha_{\beta, \gamma} = 2 \bigg( \dfrac{e \gamma}{\beta} \bigg)^{\beta / \gamma}$$
is a normalization constant, and $H(\omega)$ is the Heaviside step function.
Besides showing an additional degree of freedom in the GMWs, the two parameters $\beta$ and $\gamma$ control the time-domain and frequency-domain
decay respectively.

Example of Morse wavelet transform on the doppler test function.

=#

# Basic usage example on a doppler test function.

n=2047;
t = range(0,n/1000,length=n); # 1kHz sampling rate
f = testfunction(n, "Doppler");
p1=plot(t, f,legend=false,title="Doppler",xticks=false)
c = wavelet(Morlet(π), β=2);
res = ContinuousWavelets.cwt(f, c)
# plotting
freqs = getMeanFreq(ContinuousWavelets.computeWavelets(n, c)[1])
freqs[1] = 0
p2=heatmap(t,freqs, abs.(res)', xlabel= "time (s)", ylabel="frequency (Hz)",colorbar=false)
l=@layout [a{.3h};b{.7h}]
p = plot(p1,p2,layout=l)
savefig("morse_pic/Doppler_morlet.svg")


# Example of Morse wavelet transform on the Bumps test function.

n=2047;
t = range(0,n/1000,length=n); # 1kHz sampling rate
f = testfunction(n, "Doppler");
p1=plot(t, f,legend=false,title="Doppler",xticks=false)
plot_append = []

beta_max = 4
gamma_max = 3
for i in 1:beta_max
    for j = 1:gamma_max
        c = wavelet(Morse(i+1, j, 1));
        res = ContinuousWavelets.cwt(f, c)
        # plotting
        freqs = getMeanFreq(ContinuousWavelets.computeWavelets(n, c)[1])
        freqs[1] = 0
        p2=heatmap(t,freqs, abs.(res)',title = "(ß, γ) = (" * string(i+1) *", "* string(j) * ")", colorbar=false)
        if (j == floor((1 + gamma_max) / 2)) & (i == beta_max)
            p2=heatmap(t,freqs, abs.(res)',title = "(ß, γ) = (" * string(i+1) *", "* string(j) * ")", colorbar=false, xlabel= "time (s)")
        end
        push!(plot_append, p2)
    end
end
p = plot((plot_append[i] for i in 1:beta_max * gamma_max)..., layout = (beta_max, gamma_max), size = (900, 500))
savefig("morse_pic/Doppler_morse.svg")


# Example of Morse wavelet transform on the Bumps test function.

n=2047;
f = testfunction(n, "Bumps");
p1 = plot(f, legend = false, title = "Bumps", xlims = (0, 2000), linewidth = 2)
i = 2
j = 1
c  = wavelet(Morse(i, j, 1));
res = ContinuousWavelets.cwt(f, c)
# dropping the middle peaks
res[620:1100, :] .= 0
# and smoothing the remaining peaks
res[:, 6:end] .= 0
freqs = ContinuousWavelets.getMeanFreq(f, c)
p2 = heatmap(1:n, freqs, abs.(res)', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)
dropped = real.(ContinuousWavelets.icwt(res, c, DualFrames()))
p1 = plot(f, legend = false, title = "Smoothing & dropping bumps with " * "(ß, γ) = (" * string(i) *", "* string(j) * ")", linewidth = 2)
plot!(dropped, linewidth = 3)
l = @layout [a{0.3h}; b{0.7h}]
p = plot(p1, p2, layout = l)
savefig("morse_pic/Bumps_morse_gamma_1.svg")


n=2047;
f = testfunction(n, "Bumps");
p1 = plot(f, legend = false, title = "Bumps", xlims = (0, 2000), linewidth = 2)
i = 2
j = 2
c  = wavelet(Morse(i, j, 1));
res = ContinuousWavelets.cwt(f, c)
# dropping the middle peaks
res[620:1100, :] .= 0
# and smoothing the remaining peaks
res[:, 6:end] .= 0
freqs = ContinuousWavelets.getMeanFreq(f, c)
p2 = heatmap(1:n, freqs, abs.(res)', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)
dropped = real.(ContinuousWavelets.icwt(res, c, DualFrames()))
p1 = plot(f, legend = false, title = "Smoothing & dropping bumps with " * "(ß, γ) = (" * string(i) *", "* string(j) * ")", linewidth = 2)
plot!(dropped, linewidth = 3)
l = @layout [a{0.3h}; b{0.7h}]
p = plot(p1, p2, layout = l)
savefig("morse_pic/Bumps_morse_gamma_2.svg")


# Example of Morlet wavelet transform on the HeaviSine test function.

n = 2047
f = testfunction(n, "HeaviSine")
wave = wavelet(morl)
fCWT = ContinuousWavelets.cwt(f, wave)
fRecon = ContinuousWavelets.icwt(fCWT, wave, DualFrames())
plot([real.(fRecon) f], labels = ["reconstructed" "original"])
plot(fRecon)
deNoiseCWT = copy(fCWT)
deNoiseCWT[abs.(fCWT)/norm(fCWT, Inf).<0.07] .= 0
heatmap((abs.(fCWT) / norm(fCWT, Inf))' .< 0.01)
heatmap(abs.(fCWT'), c = :viridis)
heatmap(abs.(deNoiseCWT'), c = :viridis)
heatmap(abs.(fCWT' - deNoiseCWT'))
fRecon = ContinuousWavelets.icwt(deNoiseCWT, wave, DualFrames())
plot([fRecon f], labels = ["reconstructed" "original"])
plot(plot(f, title = "with noise"), plot(fRecon, title = "threshold and inverted"), layout = (2, 1), legend = false)
fRecon = ContinuousWavelets.icwt(fCWT, wave, DualFrames())
p = plot([real.(fRecon) f])
savefig("morse_pic/HeaviSine_morlet.svg")


# Example of Morse wavelet transform on the HeaviSine test function.

n = 2047
f = testfunction(n, "HeaviSine")
wave = wavelet(Morse(2, 3, 1));
fCWT = ContinuousWavelets.cwt(f, wave)
fRecon = ContinuousWavelets.icwt(fCWT, wave, DualFrames())
plot([real.(fRecon) f], labels = ["reconstructed" "original"])
plot(fRecon)
deNoiseCWT = copy(fCWT)
deNoiseCWT[abs.(fCWT)/norm(fCWT, Inf).<0.07] .= 0
heatmap((abs.(fCWT) / norm(fCWT, Inf))' .< 0.01)
heatmap(abs.(fCWT'), c = :viridis)
heatmap(abs.(deNoiseCWT'), c = :viridis)
heatmap(abs.(fCWT' - deNoiseCWT'))
fRecon = ContinuousWavelets.icwt(deNoiseCWT, wave, DualFrames())
plot([fRecon f], labels = ["reconstructed" "original"])
plot(plot(f, title = "with noise"), plot(fRecon, title = "threshold and inverted"), layout = (2, 1), legend = false)
fRecon = ContinuousWavelets.icwt(fCWT, wave, DualFrames())
p = plot([real.(fRecon) f])
savefig("morse_pic/HeaviSine_morse.svg")


# Example of Morlet wavelet transform on the Blocks test function.

n = 2047
f = testfunction(n, "Blocks") + 0.10randn(n)
wave = wavelet(morl)
fCWT = ContinuousWavelets.cwt(f, wave)
fRecon = ContinuousWavelets.icwt(fCWT, wave, DualFrames())
plot([real.(fRecon) f], labels = ["reconstructed" "original"])
plot(fRecon)
deNoiseCWT = copy(fCWT)
deNoiseCWT[abs.(fCWT)/norm(fCWT, Inf).<0.07] .= 0
heatmap((abs.(fCWT) / norm(fCWT, Inf))' .< 0.01)
heatmap(abs.(fCWT'), c = :viridis)
heatmap(abs.(deNoiseCWT'), c = :viridis)
heatmap(abs.(fCWT' - deNoiseCWT'))
fRecon = ContinuousWavelets.icwt(deNoiseCWT, wave, DualFrames())
plot([fRecon f], labels = ["reconstructed" "original"])
plot(plot(f, title = "with noise"), plot(fRecon, title = "threshold and inverted"), layout = (2, 1), legend = false)
fRecon = ContinuousWavelets.icwt(fCWT, wave, DualFrames())
p = plot([real.(fRecon) f], ylim = [-3, 5])
savefig("morse_pic/Blocks_morlet.svg")


# Example of Morse wavelet transform on the Blocks test function.

n = 2047
f = testfunction(n, "Blocks") + 0.10randn(n)
wave = wavelet(Morse(1.2, 2, 1));
fCWT = ContinuousWavelets.cwt(f, wave)
fRecon = ContinuousWavelets.icwt(fCWT, wave, DualFrames())
plot([real.(fRecon) f], labels = ["reconstructed" "original"])
plot(fRecon)
deNoiseCWT = copy(fCWT)
deNoiseCWT[abs.(fCWT)/norm(fCWT, Inf).<0.07] .= 0
heatmap((abs.(fCWT) / norm(fCWT, Inf))' .< 0.01)
heatmap(abs.(fCWT'), c = :viridis)
heatmap(abs.(deNoiseCWT'), c = :viridis)
heatmap(abs.(fCWT' - deNoiseCWT'))
fRecon = ContinuousWavelets.icwt(deNoiseCWT, wave, DualFrames())
plot([fRecon f], labels = ["reconstructed" "original"])
plot(plot(f, title = "with noise"), plot(fRecon, title = "threshold and inverted"), layout = (2, 1), legend = false)
fRecon = ContinuousWavelets.icwt(fCWT, wave, DualFrames())
p = plot([real.(fRecon) f], ylim = [-3, 5])
savefig("morse_pic/Blocks_morse.svg")


# Example of Morlet wavelet transform on the test function overall.

n=2047;
f = cat(testfunction(n, "Doppler"), testfunction(n, "Blocks"), testfunction(n, "Bumps"), testfunction(n, "HeaviSine"), dims = 2)

i = 3
j = 1
c  = wavelet(morl);
res = abs.(circshift(ContinuousWavelets.cwt(f, c), (0,1,0)))
freqs = ContinuousWavelets.getMeanFreq(f, c)
p1 = plot(f[:,1], legend = false, title = "Doppler", linewidth = 2)
p2 = plot(f[:,2], legend = false, title = "Blocks", linewidth = 2)
p3 = plot(f[:,3], legend = false, title = "Bumps", linewidth = 2)
p4 = plot(f[:,4], legend = false, title = "HeaviSine", linewidth = 2)

freqs = ContinuousWavelets.getMeanFreq(f, c)
h1 = heatmap(1:n, freqs, abs.(res[:,:,1])', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)
h2 = heatmap(1:n, freqs, abs.(res[:,:,2])', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)
h3 = heatmap(1:n, freqs, abs.(res[:,:,3])', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)
h4 = heatmap(1:n, freqs, abs.(res[:,:,4])', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)

p = plot(p1, p2, p3, p4, h1, h2, h3, h4, layout = (2,4), size = (800, 500))
savefig("morse_pic/test_overall_morlet.svg")


# Example of Morse wavelet transform on the test function overall.

n=2047;
f = cat(testfunction(n, "Doppler"), testfunction(n, "Blocks"), testfunction(n, "Bumps"), testfunction(n, "HeaviSine"), dims = 2)

i = 3
j = 2
c  = wavelet(Morse(i, j, 1));
res = abs.(circshift(ContinuousWavelets.cwt(f, c), (0,1,0)))
freqs = ContinuousWavelets.getMeanFreq(f, c)
p1 = plot(f[:,1], legend = false, title = "Doppler", linewidth = 2)
p2 = plot(f[:,2], legend = false, title = "Blocks", linewidth = 2)
p3 = plot(f[:,3], legend = false, title = "Bumps", linewidth = 2)
p4 = plot(f[:,4], legend = false, title = "HeaviSine", linewidth = 2)

freqs = ContinuousWavelets.getMeanFreq(f, c)
h1 = heatmap(1:n, freqs, abs.(res[:,:,1])', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)
h2 = heatmap(1:n, freqs, abs.(res[:,:,2])', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)
h3 = heatmap(1:n, freqs, abs.(res[:,:,3])', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)
h4 = heatmap(1:n, freqs, abs.(res[:,:,4])', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)

p = plot(p1, p2, p3, p4, h1, h2, h3, h4, layout = (2,4), size = (800, 500))
savefig("morse_pic/test_overall_morse.svg")



