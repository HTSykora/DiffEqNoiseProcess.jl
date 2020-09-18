@testset "ConstructedWienerGrid" begin

  using DiffEqNoiseProcess, DiffEqBase, Test

  t = 0:0.001:1
  WZWiener_trig = generate_ConstructedWienerGrid(t,Nh=10,fnspace = Trigonometric_LCConstructor)
  WZWiener_fourier = generate_ConstructedWienerGrid(t,Nh=10,fnspace = Fourier_KLConstructor)
  iterated_Strat_delay(WZWiener_trig,0.5,0.1,0.01) # W, t, τ, dt
  iterated_Strat_delay(WZWiener_fourier,0.5,0.1,0.01) # W, t, τ, dt
  double_Strat_highres(WZWiener_trig,0.5,0.01) # W, t, dt
  double_Strat_highres(WZWiener_fourier,0.5,0.01) # W, t, dt
  
  dt = 0.1
  calculate_step!(WZWiener_trig,dt,nothing,nothing)
  calculate_step!(WZWiener_fourier,dt,nothing,nothing)

  for i in 1:10
    accept_step!(WZWiener_trig,dt,nothing,nothing)
    accept_step!(WZWiener_fourier,dt,nothing,nothing)
  end

  W = generate_ConstructedWienerGrid(t, Nh=10, fnspace = Fourier_KLConstructor)
  prob = NoiseProblem(W,(0.0,1.0))
  sol = solve(prob;dt=0.1)
  iterated_Strat_delay(W,0.5,0.1,0.01) # W, t, τ, dt
  double_Strat_highres(W,0.5,0.01) # W, t, dt

  @test !sol.step_setup
  @test_throws ErrorException accept_step!(sol,dt,nothing,nothing)

  W = generate_ConstructedWienerGrid(t, W0 = zeros(8), Nh=10) # default: fnspace = Fourier_KLConstructor
  prob = NoiseProblem(W,(0.0,1.0))
  sol = solve(prob;dt=0.1)
  iterated_Strat_delay(W,0.5,0.1,0.01) # W, t, τ, dt
  double_Strat_highres(W,0.5,0.01) # W, t, dt

  dt = 1//1000
  t = 0:dt:1
  W = generate_ConstructedWienerGrid(t)
  prob_rational = NoiseProblem(W,(0,1))
  sol = solve(prob_rational; dt=1//10)
  iterated_Strat_delay(W,0.5,0.1,0.01) # W, t, τ, dt
  double_Strat_highres(W,0.5,0.01) # W, t, dt
end
