module EspressoSimulations

using LinearAlgebra
using DifferentialEquations
using NumericalIntegration: cumul_integrate
using Trapz: trapz
# using Interpolations

# ODE Setup
export EspressoParameters, SimulationParameters, to_odeproblem

# process ODESolution
export extraction_yield, extraction_percentage

include("params.jl")
include("model.jl")
include("process.jl")

end
