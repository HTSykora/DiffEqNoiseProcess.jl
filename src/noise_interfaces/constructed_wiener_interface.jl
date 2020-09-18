function save_noise!(W::ConstructedWienerGrid)

end

# Interpolation of the constructed Wiener process
function WZ_interpolant(t, t0, t1, W0::WType, WZ::WienerFnSpace{Nh,WType}) where {Nh,WType<:Number}
  W0 + WZ(t - t0, t1 - t0)
end

function WZ_interpolant(t, t0, t1, W0::WType, WZ::WienerFnSpace{Nh,WType}) where {Nh,WType<:AbstractArray}
  out = similar(W0);
  WZ_interpolant!(out,t, t0, t1, W0, WZ)
  out
end

function WZ_interpolant!(out::WType,t, t0, t1, W0::WType, WZ::WienerFnSpace{Nh,WType}) where {Nh,WType<:AbstractArray}
  WZ(out, t - t0, t1 - t0)
  out .= W0 .+ out
  out
end

function interpolate!(W::ConstructedWienerGrid, t)
  ts, timeseries, fspace = W.t, W.W, W.fspace
  sign(W.dt) * t > sign(W.dt) * ts[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Build a longer NoiseGrid to cover the integration.")
  sign(W.dt) * t < sign(W.dt) * ts[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Build a longer NoiseGrid to cover the integration.")
  tdir = sign(ts[end] - ts[1])

  if t isa Union{Rational,Integer}
      @inbounds i = searchsortedfirst(ts, t, rev=tdir < 0) # It's in the interval ts[i-1] to ts[i]
  else
      @inbounds i = searchsortedfirst(ts, t - tdir * 10eps(typeof(t)), rev=tdir < 0)
  end

  @inbounds if (t isa Union{Rational,Integer} && ts[i] == t) || (isapprox(t, ts[i]; atol=100eps(typeof(t)), rtol=100eps(t)))
      val1 = timeseries[i]
  elseif ts[i - 1] == t # Can happen if it's the first value!
      val1 = timeseries[i - 1]
  else
      dt = ts[i] - ts[i - 1]
      Θ = (t - ts[i - 1]) / dt
      val1 = WZ_interpolant(t,ts[i-1],ts[i], timeseries[i - 1], fspace[i - 1])
  end
  val1, nothing
end

function interpolate!(out1,out2,W::ConstructedWienerGrid, t)
  ts, timeseries, fspace = W.t, W.W, W.fspace
  sign(W.dt) * t > sign(W.dt) * ts[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Build a longer NoiseGrid to cover the integration.")
  sign(W.dt) * t < sign(W.dt) * ts[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Build a longer NoiseGrid to cover the integration.")
  tdir = sign(ts[end] - ts[1])

  if t isa Union{Rational,Integer}
      @inbounds i = searchsortedfirst(ts, t, rev=tdir < 0) # It's in the interval ts[i-1] to ts[i]
  else
      @inbounds i = searchsortedfirst(ts, t - tdir * 10eps(typeof(t)), rev=tdir < 0)
  end

  @inbounds if (t isa Union{Rational,Integer} && ts[i] == t) || (isapprox(t, ts[i]; atol=100eps(typeof(t)), rtol=100eps(t)))
      val1 = timeseries[i]
  elseif ts[i - 1] == t # Can happen if it's the first value!
      val1 = timeseries[i - 1]
  else
      dt = ts[i] - ts[i - 1]
      Θ = (t - ts[i - 1]) / dt
      val1 = WZ_interpolant(t,ts[i-1],ts[i], timeseries[i - 1], fspace[i - 1])
  end
  @inbounds if (t isa Union{Rational,Integer} && ts[i] == t) || (isapprox(t, ts[i]; atol = 100eps(typeof(t)), rtol = 100eps(t)))
    copyto!(out1,timeseries[i])
  elseif ts[i-1] == t # Can happen if it's the first value!
    copyto!(out1,timeseries[i-1])
  else
    WZ_interpolant!(out1,t,ts[i-1],ts[i], timeseries[i - 1], fspace[i - 1])
  end
end


function calculate_step!(W::ConstructedWienerGrid,dt,u,p)
  t = W.curt+dt
  if typeof(t) <: AbstractFloat && abs(t - W.t[end]) < 100eps(typeof(dt))
    t = W.t[end]
  end
  if isinplace(W)
    interpolate!(W.dW,W.dZ,W,t)
    W.dW .-= W.curW
  else
    new_W, _ = W(t)
    W.dW = new_W - W.curW
  end
  W.dt = dt
end

function accept_step!(W::ConstructedWienerGrid,dt,u,p,setup_next=true)
  W.step_setup == false && error("Stepped past the defined domain for the NoiseGrid")

  if isinplace(W)
    W.curW .+= W.dW
  else
    W.curW += W.dW
  end
  W.curt += W.dt

  W.dt = dt #dtpropose
  if sign(W.dt)*(W.curt + W.dt) > sign(W.dt)*W.t[end]
    setup_next = false
    W.step_setup = false
  end

  if setup_next
    calculate_step!(W,dt,u,p)
  end
  return nothing
end

function reject_step!(W::ConstructedWienerGrid,dtnew,u,p)
  calculate_step!(W,dtnew,u,p)
  return nothing
end

function setup_next_step!(W::ConstructedWienerGrid,u,p)
  calculate_step!(W,W.dt,u,p)
  return nothing
end
