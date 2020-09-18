# Kahunen-Loéve: random Fourier series
function Base.copy(fc::RandomFourierSeries{Nh,ξType})
    # TrigonometricCONS{Nh,ξType}(copy(fc.ξs)) # TODO: is a shallow copy safe?
    deepcopy(fc) # TODO: is shallow copy safe?
end

"""
RandomFourierSeries(Nh, ξprototype = 0.0)

Create a struct for the Kahunen-Loéve: type Wong-Zakai representation of a Wiener process on the interval [0,1] using random Fourier series with `2*Nh` base functions (scalable to time interval [tₙ,tₙ₊₁]).
The object is callable as:
    `(fc::RandomFourierSeries)(t_,dt)`
where `dt = tₙ₊₁ - tₙ` and `t_ = t - tₙ`.

`ξprototype <: Union{Number, AbstractVector{<:Number}}` determines the noise type (scalar or diagonal noise supported)

# Examples
```jldoctest
julia> t = 1.5; t0 = 1.; t1 = 5.;
julia> W = RandomFourierSeries(10,0.);
julia> W(t-t0, t1-t0)
0.5772330051084346
```
"""
function RandomFourierSeries(Nh::Integer, ξprototype::ξType = 0.0, a0corr = a0_correctionterm(Nh)) where {ξType <: Number}
    ξ = randn(ξType)
    ζₖ = [randn(ξType) for k in 1:Nh]
    raw_a0 = 1/(sqrt(2)*π) * sum(ζₖ[k]/k for k in 1:Nh)  - a0corr * randn(ξType) # = a₀/(2√dt)
    RandomFourierSeries(Nh,ξ,raw_a0,ζₖ,[randn(ξType) for k in 1:Nh])
end
function RandomFourierSeries(Nh::Integer, ξprototype::ξType, a0corr = a0_correctionterm(Nh)) where {ξType <: Vector{<:ξelType}} where ξelType <: Number
    ξ = randn(ξelType,size(ξprototype))
    ζₖ = [randn(ξelType, size(ξprototype)) for k in 1:Nh]
    raw_a0 = 1/(sqrt(2)*π) .* sum(ζₖ[k]/k for k in 1:Nh)  .- a0corr .* randn(ξelType, size(ξprototype)) # a₀/2√dt
    RandomFourierSeries(Nh,ξ,raw_a0,ζₖ,[randn(ξelType, size(ξprototype)) for k in 1:Nh])
end


#General functions
# Generate function space from type
getfunctionspace(constructor::Fourier_KLConstructor,Nh,W0) = RandomFourierSeries(Nh, W0, constructor.a0corr)
getfunctionspaceType(constructor::Fourier_KLConstructor,Nh,W0::T) where T = RandomFourierSeries{Nh,T,Vector{T}}

function (fc::RandomFourierSeries{Nh})(t,dt) where Nh
    sqdt = sqrt(dt);
     t / sqdt * fc.ξ .+ sqdt .* (fc.raw_a0 .+ sum((fc.ζₖ[k] * cos(t * (2*k*π/dt)) + fc.ηₖ[k] * sin(t * (2*k*π/dt))) / k for k in 1:Nh)/(sqrt(2)*π))
end 
