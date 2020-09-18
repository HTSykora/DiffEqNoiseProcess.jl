# function Iτ_Strat(W, t_,tτ_,dt_ )
    # 0.5*(W(t_+dt_) - W(t_)) * (W(tτ_+dt_) - W(tτ_))
# end

#Iterated Strationovich integral for the Levy-Ciesielski representation utilizing trigonometric functions(still needs multiplication with `dt`)
# ∫₀ᵈᵗ Wτ(t-τ) dW(t) / dt for diagonal noises
function _Itτ_Strat(W::TrigonometricCONS{Nh},Wτ::TrigonometricCONS{Nh}) where Nh
    Itτ = W.ξ₁*Wτ.ξ₁/2
    for k in 1:Nh
        Itτ += sqrt(2)*(W.ξ₁ .* Wτ.ξ₂ₖ[k] - Wτ.ξ₁ .* W.ξ₂ₖ[k])/(2k*π)
        Itτ += (W.ξ₂ₖ[k] .* Wτ.ξ₂ₖ₊₁[k] - Wτ.ξ₂ₖ[k] .* W.ξ₂ₖ₊₁[k])/(2k*π)
    end

    Itτ
end


#Iterated Strationovich integral for the random Fourier series representation (still needs multiplication with `dt`)
function _Itτ_Strat(W::RandomFourierSeries{Nh},Wτ::RandomFourierSeries{Nh}) where Nh # ∫₀ᵈᵗ Wτ(t-τ) dW(t) / dt for diagonal noises
    Itτ = (W.ξ * Wτ.ξ)/2 - (W.raw_a0 .* Wτ.ξ - Wτ.raw_a0 .* W.ξ) # TODO: check raw_a0 behavior!!!
    for k in 1:Nh
        Itτ += (1/(2k*π)) .* (Wτ.ζₖ[k] .* W.ηₖ[k] .- W.ζₖ[k] .* Wτ.ηₖ[k])
    end

    Itτ
end


function Itτ_Strat(W::ConstructedWienerGrid{T,N,Tt,T2,<:AbstractVector{<:WienerFnSpace{<:Nh}}}, t_,tτ_,dt_ ) where {T,N,Tt,T2,Nh}
    @assert t_ >= tτ_ "The condition t >= t-τ is not satisfied!"
    # if (isapprox(t_, tτ_; atol=100eps(typeof(t_)), rtol=100eps(t_)))
    #     return W.dW^2
    # end

    ts, timeseries, fspace = W.t, W.W, W.fspace
    tdir = sign(ts[end] - ts[1])
    @assert tdir > 0 "Method for backwards integration is not impelemented!"

    sign(dt_) * t_ > sign(dt_) * ts[end] || sign(dt_) * tτ_+dt_ > sign(dt_) * ts[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Build a longer NoiseGrid to cover the integration.")
    if sign(dt_) * t_ <= sign(dt_) * ts[1]
        return zero(W.W[1])
    end
    if sign(dt_) * tτ_ <= sign(dt_) * ts[1] && sign(dt_) * tτ_ + dt_ <= sign(dt_) * ts[1]
        return zero(W.W[1])
    elseif sign(dt_) * tτ_ <= sign(dt_) * ts[1] && sign(dt_) * tτ_ + dt_ > sign(dt_) * ts[1]
        tτ = ts[1]
        t = t_ + (ts[1] - tτ_)
        dt = dt_ - (ts[1] - tτ_)
    else
        t = t_
        tτ = tτ_
        dt = dt_
    end

    i0,bi0 = gettimeidx(t,ts,tdir)
    i1,bi1 = gettimeidx(t+dt,ts,tdir)
    i1 -= 1
    j0,bj0 = gettimeidx(tτ,ts,tdir)
    j1,bj1 = gettimeidx(tτ+dt,ts,tdir)
    j1 -= 1

    if !any([bi0,bi1,bj0,bj1])
        if ts[i0 - 1] == t # Can happen if it's the first value!
            i0 = i0-1
        elseif ts[j0 - 1] == t # Can happen if it's the first value!
            j0 = j0-1
        else
            error("No interpolation is possible in case of iterated stochastic integrals")
        end
    end
    @assert j1 - j0 == i1 - i0 "t = ts[i] or t-τ = ts[j] is wrongly selected"
    Δij = i0 - j0;
    dt * sum(_Itτ_Strat(fspace[i],fspace[i-Δij]) for i in i0:i1) / sqrt(i1-i0+1) # TODO: Multiple timesteps??? WHY SQRT???
end
  
function I_Strat(W::ConstructedWienerGrid{T,N,Tt,T2,<:AbstractVector{<:WienerFnSpace{<:Nh}}}, t_,dt_ ) where {T,N,Tt,T2,Nh}

    ts, timeseries, fspace = W.t, W.W, W.fspace
    tdir = sign(ts[end] - ts[1])
    @assert tdir > 0 "Method for backwards integration is not impelemented!"

    sign(dt_) * t_ > sign(dt_) * ts[end]  && error("Solution interpolation cannot extrapolate past the final timepoint. Build a longer NoiseGrid to cover the integration.")
    if sign(dt_) * t_ < sign(dt_) * ts[1]
        error("Solution interpolation cannot extrapolate before the first timepoint. Build a longer NoiseGrid to cover the integration.")
    end
    t = t_
    dt = dt_

    i0,bi0 = gettimeidx(t,ts,tdir)
    i1,bi1 = gettimeidx(t+dt,ts,tdir)
    i1 -= 1

    if !any([bi0,bi1])
        if ts[i0 - 1] == t # Can happen if it's the first value!
            i0 = i0-1
        else
            error("No interpolation is possible in case of iterated stochastic integrals")
        end
    end
    dt * sum((getincrementvariable(fspace[i]) .^ 2)*0.5 for i in i0:i1) / sqrt(i1-i0+1) # TODO: Multiple timesteps???
end

function gettimeidx(t,ts,tdir)
    if t isa Union{Rational,Integer}
        @inbounds i = searchsortedfirst(ts, t, rev=tdir < 0)  
    else
        @inbounds i = searchsortedfirst(ts, t - tdir * 10eps(typeof(t)), rev=tdir < 0)
    end
    @inbounds  bi = (t isa Union{Rational,Integer} && ts[i] == t) || (isapprox(t, ts[i]; atol=100eps(typeof(t)), rtol=100eps(t)))
    i,bi
end