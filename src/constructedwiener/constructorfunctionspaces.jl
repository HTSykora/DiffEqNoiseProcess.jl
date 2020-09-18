# Complete OrthoNormal Systems (CONS) used for the Wong-Zakai construction of W(t) 
# Candidates:
# Levy-Ciesielski: Trigonometric, Haar, Shauder
# Kahunen-Loéve: random Fourier series
 
abstract type WienerConstructor{Nh} end

# Levy-Ciesielski
"""
Trigonometric_LCConstructor(Nh)

Type for the Levy-Ciesielski type Wong-Zakai representation of a Wiener process on the interval [tₙ,tₙ₊₁] using a linear and `Nh` trigonometric functions.
"""
struct Trigonometric_LCConstructor{Nh} <: WienerConstructor{Nh} end
Trigonometric_LCConstructor(Nh::Integer) = Trigonometric_LCConstructor{Nh}()

# Nh: cardinality of the CONS used 
# ξ: stochastic coefficients of the decomposition
abstract type WienerFnSpace{Nh,ξType} end 

struct TrigonometricCONS{Nh,ξType,ξVType} <: WienerFnSpace{Nh,ξType} # Trigonomietric complete orthonormal system
    ξ₁::ξType #  ξconst
    ξ₂ₖ::ξVType # ξ_{2k} 
    ξ₂ₖ₊₁::ξVType # ξ_{2k+1} 
    TrigonometricCONS(Nh,ξ₁::ξType,ξ₂ₖ::ξVType,ξ₂ₖ₊₁::ξVType) where ξVType<:AbstractVector{<:ξType} where ξType = new{Nh,typeof(ξ₁), typeof(ξ₂ₖ)}(ξ₁,ξ₂ₖ,ξ₂ₖ₊₁)
end
getincrementvariable(fnspace::TrigonometricCONS) = fnspace.ξ₁

# Kahunen-Loéve:
a0_correctionterm(Nh) = sqrt(1/12 - (1/(2*(π^2)))*sum(1/k^2 for k in 1:Nh))

"""
Fourier_KLConstructor(Nh)

Type for the Kahunen-Loéve type Wong-Zakai representation of a Wiener process on the interval [tₙ,tₙ₊₁] using random Fourier series with `2*Nh` base functions.
"""
struct Fourier_KLConstructor{Nh} <: WienerConstructor{Nh}
    a0corr::Float64
end

Fourier_KLConstructor(Nh::Integer) = Fourier_KLConstructor{Nh}(a0_correctionterm(Nh))

# Origin: random Fourier series of the Brownian Bridge B(t) = W(t) - t/Δ W(Δ) t∈[0,Δ]
struct RandomFourierSeries{Nh,Type,VType} <: WienerFnSpace{Nh,Type}
    # a0corr::Type # a0 correction term for the simulations/integrals?
    ξ::Type
    raw_a0::Type # = a₀/(2√dt)
    ζₖ::VType # coefficient of aₖ in aₖ cos(2kπ/dt*t)
    ηₖ::VType # generator of bₖ in bₖ sin(2kπ/dt*t)
    
    RandomFourierSeries(Nh,ξ::Type,raw_a0::Type,ζₖ::VType,ηₖ::VType) where VType<:AbstractVector{<:Type} where Type = new{Nh,typeof(raw_a0), typeof(ζₖ)}(ξ,raw_a0,ζₖ,ηₖ)
end
getincrementvariable(fnspace::RandomFourierSeries) = fnspace.ξ


# GenerateFunctions

function generate_ConstructedWienerGrid(t; W0=0.0, Nh=4, fnspace = Fourier_KLConstructor)
    # t = t0:dt(:T + 100eps(T));
    
    W, CONSs = generate_ConstructedWienerGrid(W0, fnspace(Nh), t[2] - t[1], t)
    ConstructedWienerGrid(t, W, CONSs)
end

function generate_ConstructedWienerGrid(W0::WType, fnspace::fnType, dt, t) where {WType,fnType <: WienerConstructor{<:Nh}} where {Nh}
    W = Vector{WType}(undef, 0); sizehint!(W, length(t))
    fnspaceV = Vector{getfunctionspaceType(fnspace, Nh, W0)}(undef, 0); sizehint!(fnspaceV, length(t) - 1);
    sqdt = sqrt(dt);

    push!(W, W0);
    push!(fnspaceV, getfunctionspace(fnspace, Nh, W0));
    push!(W, W0 .+ sqdt * getincrementvariable(fnspaceV[end])) # dt = t[2]-t[1] ?
    
    for enumerate(i, t) in t[3:end]
        push!(fnspaceV, getfunctionspace(fnspace, Nh, W0));
        push!(W, W[end] .+ sqdt .* getincrementvariable(fnspaceV[end])) # dt = t-t[i+1] ?
    end
    W, fnspaceV
end