"""
    nOctaves, totalWavelets, sRanges, sWidth = getNWavelets(n1,c)

utility for understanding the spacing of the wavelets. `sRanges` is a list of
the s values used in each octave. sWidth is a list of the corresponding
variance adjustments
"""
function getNWavelets(n1, c::CWT)
    nOctaves = getNOctaves(n1, c)
    n,nSpace = setn(n1,c)
    isAve = !(typeof(c.averagingType) <: NoAve)
    if nOctaves ≤ c.averagingLength + getMinScaling(c) # there isn't enough space for anything but averaging.
        totalWavelets = isAve
        sRanges = [1.0]
        sWidth = [1.0]
        return nOctaves, totalWavelets, sRanges, sWidth
    end
    sRange = 2 .^ (polySpacing(nOctaves, c))
    totalWavelets = round(Int, length(sRange) + isAve)
    sWidth = varianceAdjust(c,totalWavelets,nOctaves)
    return nOctaves, totalWavelets, sRange, sWidth
end




"""
different wavelet familes need to end at a different number of octaves because they have different
"""
getNOctaves(n1,c::CWT{W,T, M, N}) where {W, T, N, M} = log2(n1>>1+1) + c.extraOctaves
# choose the number of octaves so the last mean, which is at s*σ[1]
# is 3 standard devations away from the end
getNOctaves(n1,c::CWT{W,T, Morlet, N}) where {W, T, N} = log2((n1>>1+1)/(c.σ[1]+3)) + c.extraOctaves
getNOctaves(n1,c::CWT{W,T, <:Paul, N}) where {W, T, N} = log2((n1>>1+1)/(2c.α+5)) + c.extraOctaves
# choose the number of octaves so the last mean is 2 standard deviations from the end
function getNOctaves(n1,c::CWT{W,T, <:Dog, N}) where {W, T, N}
    μ = getMean(c)
    σ = getStd(c)
    log2(n1>>1/(μ+5σ)) + c.extraOctaves
end
# choose the number of octaves so the smallest support is twice the qmf
getNOctaves(n1,c::CWT{W,T, <:ContOrtho, N}) where {W, T, N} = log2(n1) - 2 - log2(length(qmf(c.waveType))) + c.extraOctaves

"""
As with the last octave, different wavelet families have different space decay rates, and in the case of symmetric or zero padding we don't want wavelets that bleed across the boundary from the center.
"""
getMinScaling(c::CWT{W,T,M,N}) where {W,T,N,M} = 0 # by default all scales are allowed (all of the orthogonal transforms)
getMinScaling(c::CWT{W,T,<:Morlet,N}) where {W,T,N} = 1/(c.β)^.8 # morlet is slightly too large at the boundaries by default
getMinScaling(c::CWT{W,T,<:Paul,N}) where {W,T,N} = 2/(2c.α+1) # Paul presents some difficulties, as the decay changes quickly (like 1/t^(α+1))
getMinScaling(c::CWT{W,T,<:Dog,N}) where {W,T,N} = 2 # like morlet, the decay for Dog is exponential and consistent across derivatives



function varianceAdjust(this::CWT{W,T, M, N}, totalWavelets, nOct) where {W,T,N, M}
    # increases the width of the wavelets by σ[i] = (1+a(total-i)ᴾ)σₛ
    # which is a polynomial of order p going from 1 at the highest frequency to
    # sqrt(p) at the lowest
    β = this.β
    x = exp(nOct)/(100+exp(nOct)) # the fewer octaves there are, the smaller the width adjustment we need
    a = (β^x-1)/(totalWavelets-1)^β
    sWidth = 1 .+a .*(totalWavelets .- (1:totalWavelets)).^β
    if any(isnan.(sWidth))
        return [1]
    end
    return sWidth
end


function polySpacing(nOct, c)
    a = getMinScaling(c) + c.averagingLength
    O = nOct
    β = c.β
    Q = c.Q
    if O ≤ a
        # averagingLength and the min scaling are too large for anything to be done
        return [1.0]
    elseif 0 < O - a ≤ 1
        # there's only one octave, just return the linear rate of Q
        return range(a, O, length=1+round(Int, Q))
    end
    # x is the index, y is the scale
    # y= b*x^(1/β), solve for a and b with
    # (x₀,y₀)=(1, aveLength + minScale)
    # dy/dx = 1/s, so the Quality factor gives the slope at the last frequency
    b = (β/Q)^(1/β) * (O)^((β-1)/β)
    # the point x so that the second condition holds
    lastWavelet = Q * (O)/β
    # the point so that the first wavelet is at a
    firstWavelet = (a/b)^β
    # step size so that there are actually Q wavelets in the last octave
    startOfLastOctave = ((nOct-1)/b)^β
    stepSize = (lastWavelet - startOfLastOctave)/Q
    samplePoints = range(firstWavelet, stop=lastWavelet,
                         step = stepSize)#, length = round(Int, O * Q^(1/p)))
    return b .* (samplePoints).^(1/β)
end

# a utility to just get the start, stop, and step size used in polySpacing. Only used for explanatory purposes
function genSamplePoints(nOct, c)
    a = getMinScaling(c) + c.averagingLength
    O = nOct
    β = c.β
    Q = c.Q
    if O ≤ a
        # averagingLength and the min scaling are too large for anything to be done
        return [1.0]
    elseif 0 < O - a ≤ 1
        # there's only one octave, just return the linear rate of Q
        return range(a, O, length=1+round(Int, Q))
    end
    # x is the index, y is the scale
    # y= b*x^(1/β), solve for a and b with
    # (x₀,y₀)=(1, aveLength + minScale)
    # dy/dx = 1/s, so the Quality factor gives the slope at the last frequency
    b = (β/Q)^(1/β) * (O)^((β-1)/β)
    # the point x so that the second condition holds
    lastWavelet = Q * (O)/β
    # the point so that the first wavelet is at a
    firstWavelet = (a/b)^β
    # step size so that there are actually Q wavelets in the last octave
    startOfLastOctave = ((nOct-1)/b)^β
    stepSize = (lastWavelet - startOfLastOctave)/Q
    firstWavelet, lastWavelet, stepSize
end

"""

adjust the length of the storage based on the boundary conditions
"""
function setn(n1, c)
    if boundaryType(c) <: ZPBoundary
        base2 = ceil(Int,log2(n1 + 1));   # power of 2 nearest to n1
        nSpace = 2^(base2)
        n = nSpace>>1 + 1
    elseif boundaryType(c) <: SymBoundary
        # n1+1 rather than just n1 because this is going to be used in an rfft
        # for real data
        n = n1 + 1
        nSpace = 2*n1
    elseif boundaryType(c) <: PerBoundary
        n= n1>>1 + 1
        nSpace = n1
    else
        error("No such boundary as $(boundaryType(c))")
    end
    return n, nSpace
end






# computing averaging function utils

# this makes sure the averaging function is real and positive
adjust(c::CWT) = 1
adjust(c::CWT{W, T, <:Dog}) where {W,T} = 1/(im)^(c.α)

function locationShift(c::CWT{W, T, <:Morlet, N}, s, ω, sWidth) where {W,T,N}
        s0 = 3*c.σ[1]/4 *s*sWidth
        ω_shift = ω .+ c.σ[1] * s0
    return (s0, ω_shift)
end

function locationShift(c::CWT{W, T, <:Dog, N}, s, ω, sWidth) where {W,T,N}
    μNoS = getMean(c)
    μLast = μNoS * (2s)
    s0 = μLast / 2 / getStd(c)
    μ = μNoS * s0
    ω_shift = ω .+ μ
    return (s0, ω_shift)
end

"""
get the mean of the mother wavelet where s=1
"""
function getMean(c::CWT{W, T, <:Dog},s=1) where {W,T}
    Gm1 = gamma((c.α+1)/2)
    Gm2 = gamma((c.α+2)/2)
    Gm2 / Gm1 * sqrt(2) * s^2
end
getMean(c::CWT{W,T,<:Paul}, s=1) where {W,T} = (c.α+1)*s
function getMean(c::CWT{W, T, <:Morlet},s=1) where {W,T}
    return s*c.σ[1]
end
"""
    getStd(c::CWT{W, T, <:Dog}, s=1) where {W,T}
get the standard deviation of the mother wavelet
"""
getStd(c::CWT{W, T, <:Dog}, s=1) where {W,T} = sqrt(c.α + 1 - getMean(c)^2)*s^(3/2)
getStd(c::CWT{W, T, <:Paul},s=1) where {W, T} = sqrt((c.α+2)*(c.α+1)) * s
getStd(c::CWT{W, T, <:Morlet},s=1) where {W,T} = (s * c.β^.8)/sqrt(2)

function locationShift(c::CWT{W, T, <:Paul, N}, s, ω, sWidth) where {W,T,N}
    s0 = s*sqrt(c.α+1)
    ω_shift = ω .+ (c.α .+ 1) * s0
    return (s0, ω_shift)
end


function getUpperBound(c::CWT{W, T, <:Morlet, N}, s) where {W,T,N}
    return c.σ[1] * s
end

function getUpperBound(c::CWT{W, T, <:Dog, N}, s) where {W,T,N}
    return sqrt(2)* s * gamma((c.α+2)/2) / gamma((c.α+1)/2)
end

function getUpperBound(c::CWT{W, T, <:Paul, N}, s) where {W,T,N}
    return (c.α + 1) * s
end

"""
    arrayOfFreqs = getMeanFreq(Ŵ, δt=1000)
Calculate each of the mean frequencies of a collection of analytic or real wavelets Ŵ. assumes a sampling rate of 2kHz, so the maximum frequency is 1kHz. Change δt to adjust to your problem.
"""
function getMeanFreq(Ŵ, δt=2000)
    eachNorm = [norm(w,1) for w in eachslice(Ŵ,dims=ndims(Ŵ))]'
    freqs = range(0,δt/2, length= size(Ŵ,1))
    dropdims(sum(freqs .* abs.(Ŵ), dims=1) ./ eachNorm,dims=1)
end


# orthogonal wavelet utils
function getContWaveFromOrtho(c,N)
    depth, ψLen = calcDepth(c,N)
    φ, ψ, t = makewavelet(c.waveType.o, depth)
    return φ, ψ, ψLen
end

"""
    nIters, sigLength = calcDepth(w,N)
given a `CWT` with an orthogonal wavelet type, calculate the depth `nIters`
necessary to get a resulting wavelet longer than `N`. `sigLength` is the
resulting length.
"""
function calcDepth(w, N)
    origLen = length(qmf(w.waveType))
    netLen = origLen
    nIters = 0
    if N < netLen
        return nIters
    end
    while N > netLen - 2^nIters + 1
        netLen = nextLength(netLen, origLen)
        nIters += 1
    end
    return (nIters, netLen - 2^nIters + 1)
end
# the actual length calculator, as this is recursive
nextLength(curr,orig) = curr*2 + orig - 1 # upsample then convolve
# pad a vector `v` to have length `N`
padTo(v, N) = cat(v, zeros(eltype(v), max(N-length(v),0)), dims = 1)

# create interpolater for the orthogonal cases
genInterp(ψ) = interpolate(ψ, BSpline(Quadratic(Reflect(OnGrid()))))
