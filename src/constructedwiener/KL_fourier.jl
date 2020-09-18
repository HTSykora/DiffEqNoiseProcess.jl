# Kahunen-Loéve: random Fourier sereis
function Base.copy(fc::RandomFourierSeries)
    # TrigonometricCONS{Nh,ξType}(copy(fc.ξs)) # TODO: is a shallow copy safe?
    deepcopy(fc)
end

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
