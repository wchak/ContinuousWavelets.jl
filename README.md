# ContinuousWavelets

[![Build Status](https://travis-ci.com/dsweber2/ContinuousWavelets.jl.svg?branch=master)](https://travis-ci.com/dsweber2/ContinuousWavelets.jl)
[![Coverage](https://codecov.io/gh/dsweber2/ContinuousWavelets.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dsweber2/ContinuousWavelets.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://dsweber2.github.io/ContinuousWavelets.jl/dev/)

This package is an offshoot of [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl) for the continuous wavelets.
Thanks to [Felix Gerick](https://github.com/fgerick) for the initial implementation there, with extension and further adaptation by David Weber and any other contributors listed on the right.
Currently, it implements 1D continuous wavelet transforms with the following mother wavelets:

![Mothers](docs/mothers.svg)

Which covers several standard continuous wavelet families, both real and analytic, as well as continuous versions of the orthogonal wavelet transforms implemented in [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl).

## Basic Usage

Install via the package manager and load with `using`

```julia
julia> Pkg.add("ContinuousWavelets")
julia> using ContinuousWavelets
```

Basic usage example on a doppler test function.

```julia
julia> using ContinuousWavelets, Plots, Wavelets

julia> n = 2047;

julia> t = range(0, n / 1000, length=n); # 1kHz sampling rate

julia> f = testfunction(n, "Doppler");

julia> p1 = plot(t, f, legend=false, title="Doppler", xticks=false)
Plot{Plots.PyPlotBackend() n=1}

julia> c = wavelet(Morlet(π), β=2)

julia> res = ContinuousWavelets.cwt(f, c)
┌ Warning: the lowest frequency wavelet has more than 1% its max at zero, so it may not be analytic. Think carefully
│   lowAprxAnalyt = 0.06186323501016359
└ @ ContinuousWavelets ~/work/ContinuousWavelets.jl/ContinuousWavelets.jl/src/sanityChecks.jl:6
2047×31 Matrix{ComplexF64}:
 -1.48637e-6+3.8241e-19im   …  0.000109978+9.67834e-5im
 -1.48602e-6+5.15534e-19im     -8.24922e-5+0.000130656im
            ⋮               ⋱             ⋮
 0.000435175+2.30636e-19im  …  -2.47195e-6-1.97048e-8im
 0.000435027-8.28725e-19im     -2.63499e-6+4.62331e-8im
```

And now we make a scalogram to actually visualize these entries:

```julia
freqs = getMeanFreq(ContinuousWavelets.computeWavelets(n, c)[1])
freqs[1] = 0
p2 = heatmap(t, freqs, log.(abs.(res).^2)', xlabel= "time (s)", ylabel="frequency (Hz)", colorbar=false, c=cgrad(:viridis, scale=:log10))
l = @layout [a{.3h};b{.7h}]
plot(p1,p2,layout=l)
```

![Doppler](/docs/doppler.svg)

As the cwt frame is redundant, there are many choices of dual/inverse frames. There are three available in this package, `NaiveDelta()`, `PenroseDelta()`, and `DualFrames()`. As a toy example, lets knock out the middle time of the bumps function and apply a high pass filter:

```julia
f = testfunction(n, "Bumps");
p1 = plot(f, legend = false, title = "Bumps", xlims = (0, 2000), linewidth = 2)
c = wavelet(dog2, β = 2)
res = ContinuousWavelets.cwt(f, c);
# dropping the middle peaks
res[620:1100, :] .= 0
# and smoothing the remaining peaks
res[:, 10:end] .= 0
freqs = ContinuousWavelets.getMeanFreq(length(f), c)
p2 = heatmap(1:n, freqs, abs.(res)', xlabel = "time (ms)", ylabel = "Frequency (Hz)", colorbar = false, c = :viridis)
dropped = ContinuousWavelets.icwt(res, c, DualFrames());
p1 = plot(f, legend=false, title="Smoothing and dropping bumps", linewidth=2)
plot!(dropped, linewidth=3)
l = @layout [a{0.3h}; b{0.7h}]
plot(p1, p2, layout=l)
```

![Bumps](/docs/bumps.svg)

It can also handle collections of examples at the same time, should you need to do a batch of transforms:

```julia
julia> exs = cat(testfunction(n, "Doppler"), testfunction(n, "Blocks"), testfunction(n, "Bumps"), testfunction(n, "HeaviSine"), dims=2);

julia> c = wavelet(cDb2, β=2, extraOctaves=-0)


julia> res = circshift(ContinuousWavelets.cwt(exs, c), (0, 1, 0))
┌ Warning: the highest frequency wavelet has more than 1% its max at the end, so it may not be analytic. Think carefully
│   highAprxAnalyt = 0.2677814440444114
└ @ ContinuousWavelets ~/work/ContinuousWavelets.jl/ContinuousWavelets.jl/src/sanityChecks.jl:10
2047×32×4 Array{Float64, 3}:
[:, :, 1] =
 1.89367e-5  0.000266033  …  4.6727e-5    2.99983e-6
 8.33321e-5  0.000266913     1.56557e-5  -4.46419e-5
 ⋮                        ⋱  ⋮
 2.24677e-6  0.00198709   …  4.24042e-6   3.80685e-6
 2.63848e-6  0.00198004      4.3791e-6    3.47575e-6

[:, :, 2] =
  7.81007e-18  0.0226754  0.00955729  …   3.68809e-18
 -3.47114e-18  0.022684   0.00950443     -3.47114e-18
  ⋮                                   ⋱
 -9.29595e-18  0.0341512  0.0108039   …  -3.84208e-19
  1.27592e-18  0.0342157  0.0107729      -1.7043e-18

[:, :, 3] =
 -4.2736e-7   0.0059687   …  4.47839e-8  1.86209e-8
 -4.39691e-7  0.00596762     3.30771e-8  7.78201e-9
  ⋮                       ⋱  ⋮
 -9.41123e-8  0.00339924  …  8.01012e-9  4.78652e-9
 -9.36079e-8  0.0034061      8.3188e-9   4.24252e-9

[:, :, 4] =
  0.000307454  -0.0150898   -0.00391724  …   0.000301757
  6.05948e-5   -0.0152536   -0.00405883      8.45503e-5
  ⋮                                      ⋱
 -0.000307094  -0.00755439  -0.00156729  …  -0.000594673
 -0.000378125  -0.00746687  -0.00146262     -0.00051676
```

And the plot of these:

```julia
p1 = plot(plot(exs[:, 1], legend=false, title="Doppler", yticks=[], xticks=[], linewidth=2), plot(exs[:, 2], legend=false, title="Blocks", yticks=[], xticks=[], linewidth=2), plot(exs[:, 3], legend=false, title="Bumps", yticks=[], xticks=[], linewidth=2), plot(exs[:, 4], legend=false, title="HeaviSine", yticks=[], xticks=[], linewidth=2), layout=(1, 4))
p2 = plot(heatmap(identity.(res[:, :, 1])', xticks=false, yticks=[], c=:viridis, colorbar=false), heatmap(identity.(res[:, :, 2])', xticks=false, yticks=[], c=:viridis, colorbar=false), heatmap(identity.(res[:, :, 3])', xticks=false, yticks=[], c=:viridis, colorbar=false), heatmap(identity.(res[:, :, 4])', xticks=false, yticks=[], c=:viridis, colorbar=false), layout=(1, 4))
l = @layout [a{0.3h}; b{0.7h}]
plot(p1, p2, layout=l)
```

![parallel transforms](/docs/multiEx.svg)

Example of morse wavelet transform

The Generalized Morse Wavelets (GMWs) is a superfamily of analytic wavelets.

The GMWs in the frequency domain can be defined by

```math
$$\Psi_{\beta, \gamma}(\omega) := \int \psi_{\beta, \gamma}(t) e^{-i \omega t} d{t} = H(\omega) \alpha_{\beta, \gamma} \omega^\beta e^{-\omega^\gamma} $$
```

where 
```math
$\beta>0, \gamma \geq 1$ are two main parameters controlling the form of wavelet, $$\alpha_{\beta, \gamma} = 2 \bigg( \dfrac{e \gamma}{\beta} \bigg)^{\beta / \gamma}$
```
is a normalization constant, and $H(\omega)$ is the Heaviside step function.
Besides showing an additional degree of freedom in the GMWs, the two parameters $\beta$ and $\gamma$ control the time-domain and frequency-domain decay respectively.


Example of Morse wavelet transform on the doppler test function.
```julia
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
```

![Dop_morse](/docs/morse_pic/Doppler_morse.svg)

Example of Morlet wavelet transform on the HeaviSine test function.
![HeaviSine_morlet](/docs/morse_pic/HeaviSine_morlet.svg")

Example of Morse wavelet transform on the HeaviSine test function.
![HeaviSine_morse](/docs/morse_pic/HeaviSine_morse.svg")

Example of Morlet wavelet transform on the Blocks test function.
![Blocks_morlet](/docs/morse_pic/Blocks_morlet.svg")

Example of Morse wavelet transform on the Blocks test function.
![Blocks_morse](/docs/morse_pic/Blocks_morse.svg")



There are also several boundary conditions, depending on the kind of data given; the default `SymBoundary()` symmetrizes the data, while `PerBoundary()` assumes it is periodic, and `ZPBoundary` pads with zeros.
All wavelets are stored in the Fourier domain, and all transforms consist of performing an fft (possibly an rfft if the data is real) of the input, pointwise multiplication (equivalent to convolution in the time domain), and then returning to the time domain.

Perhaps somewhat unusually, the averaging function, or father wavelet, is included as an option (the bottom row for the figure above). This can be either the paired averaging function or uniform in frequency (the `Dirac` averaging). The frequency coverage of the wavelets can be adjusted both in total frequency range both below by the `averagingLength` or above by the `extraOctaves` (caveat emptor with how well they will be defined in that case). The frequency density can be adjusted both in terms of the quality/scale factor `Q`, as well as how quickly this density falls off as the frequency goes to zero via `β`. Finally, depending on what kind of norm you want to preserve, `p` determines the norm preserved in the frequency domain (so `p=1` maintains the 1-norm in frequency, while `p=Inf` maintains the 1-norm in time).

## Possible extensions

- Higher dimensional wavelets have yet to be implemented.
- A DCT implementation of the symmetric boundary to halve the space and computational costs.
