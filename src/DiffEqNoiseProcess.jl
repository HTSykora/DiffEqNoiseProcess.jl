module DiffEqNoiseProcess

using ResettableStacks, DiffEqBase, RecipesBase
using RecursiveArrayTools, StaticArrays, Random, Statistics
using LinearAlgebra, Requires

import RandomNumbers: Xorshifts

import DiffEqBase: isinplace, solve, AbstractNoiseProcess,
       DEIntegrator, AbstractNoiseProblem

import PoissonRandom, Distributions

using DiffEqBase: @..

include("init.jl")
include("types.jl")
include("wiener.jl")
include("solve.jl")
include("geometric_bm.jl")
include("compoundpoisson.jl")
include("ornstein_uhlenbeck.jl")
include("rswm.jl")
include("bridges.jl")
include("constructedwiener/constructorfunctionspaces.jl")
include("constructedwiener/LC_trigonometric.jl")
include("constructedwiener/KL_fourier.jl")
include("constructedwiener/iteratedstratonovichintegral.jl")
include("noise_interfaces/simple_noise_process_interface.jl")
include("noise_interfaces/noise_process_interface.jl")
include("noise_interfaces/noise_function_interface.jl")
include("noise_interfaces/noise_grid_interface.jl")
include("noise_interfaces/noise_approximation_interface.jl")
include("noise_interfaces/noise_wrapper_interface.jl")
include("noise_interfaces/common.jl")
include("noise_interfaces/constructed_wiener_interface.jl")
include("correlated_noisefunc.jl")

export RSWM

export NoiseProcess, SimpleNoiseProcess, adaptive_alg

export WienerProcess, WienerProcess!, SimpleWienerProcess, SimpleWienerProcess!

export RealWienerProcess, RealWienerProcess!

export BrownianBridge, BrownianBridge!

export GeometricBrownianMotionProcess, GeometricBrownianMotionProcess!

export CompoundPoissonProcess, CompoundPoissonProcess!

export GeometricBrownianBridge, GeometricBrownianBridge!

export CompoundPoissonBridge, CompoundPoissonBridge!

export OrnsteinUhlenbeckProcess, OrnsteinUhlenbeckProcess!

export NoiseWrapper, NoiseFunction, NoiseGrid, NoiseApproximation, ConstructedWienerGrid

export accept_step!, reject_step!, calculate_step!, setup_next_step!, save_noise!

export CorrelatedWienerProcess, CorrelatedWienerProcess!

export ConstructedWienerGrid, generate_ConstructedWienerGrid, Trigonometric_LCConstructor, Fourier_KLConstructor

end # module
