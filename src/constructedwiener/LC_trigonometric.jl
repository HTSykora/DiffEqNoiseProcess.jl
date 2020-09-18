# Levy-Ciesielski: Trigonometric
function Base.copy(fc::TrigonometricCONS{Nh,ξType}) where {Nh,ξType}
    # TrigonometricCONS{Nh,ξType}(copy(fc.ξs)) # TODO: is a shallow copy safe?
    deepcopy(fc)
end

function TrigonometricCONS(Nh::Integer, ξprototype::ξType = 0.0) where {ξType <: Number}
    Nh2 = Nh÷2
    TrigonometricCONS(Nh2,randn(ξType),[randn(ξType) for i in 1:Nh2],[randn(ξType) for i in 1:Nh2])
end
function TrigonometricCONS(Nh::Integer, ξprototype::ξType) where {ξType <: Vector{<:ξelType}} where ξelType <: Number
    Nh2 = Nh÷2
    TrigonometricCONS(Nh2,randn(ξelType, size(ξprototype)),
        [randn(ξelType, size(ξprototype)) for i in 1:Nh2],
        [randn(ξelType, size(ξprototype)) for i in 1:Nh2])
end


#General functions
# Generate function space from type
getfunctionspace(constructor::Trigonometric_LCConstructor,Nh,W0) = TrigonometricCONS(Nh, W0)
getfunctionspaceType(constructor::Trigonometric_LCConstructor,Nh,W0::ξT) where ξT = TrigonometricCONS{Nh÷2,ξT,Vector{ξT}}

function (fc::TrigonometricCONS{Nh})(t,dt) where Nh
    sqdt = sqrt(dt);
    t/sqdt .* fc.ξ₁ .+ sqrt(2)*sqdt/(2π) .* sum(fc.ξ₂ₖ[k]*(1-cos(2π*k/dt*t))/k + fc.ξ₂ₖ₊₁[k]*sin(2π*k/dt*t)/k for k in 1:Nh)
end 

# function generate_ConstructedWienerGrid(t; W0=0.0, Nh=4, fnspace = Trigonometric_LCConstructor)
#     # t = t0:dt:T + 100eps(T);
    
#     W, CONSs = generate_ConstructedWienerGrid(W0, fnspace(Nh), t[2] - t[1], t)
#     ConstructedWienerGrid(t, W, CONSs)
# end

# function generate_ConstructedWienerGrid(W0::WType, WZAppr::wzaType, dt, t) where {WType,wzaType <: WienerConstructor{<:Nh}} where {Nh}
#     W = Vector{WType}(undef, 0); sizehint!(W, length(t))
#     fnspaceV = Vector{getfunctionspaceType(WZAppr, Nh, W0)}(undef, 0); sizehint!(fnspaceV, length(t) - 1);
#     sqdt = sqrt(dt);

#     push!(W, W0);
#     push!(fnspaceV, getfunctionspace(WZAppr, Nh, W0));
#     push!(W, W0 .+ sqdt * fnspaceV[end].ξ₁) # dt = t[2]-t[1] ?
    
#     for enumerate(i, t) in t[3:end]
#         push!(fnspaceV, getfunctionspace(WZAppr, Nh, W0));
#         push!(W, W[end] .+ sqdt .* fnspaceV[end].ξ₁) # dt = t-t[i+1] ?
#     end
#     W, fnspaceV
# end

