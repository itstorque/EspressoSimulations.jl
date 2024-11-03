using DifferentialEquations
using NumericalIntegration
using Trapz

function extraction_percentage(sol:: ODESolution, espresso_params:: EspressoParameters, sim_params:: SimulationParameters; cumul::Bool=false)
    dimensionless_params = generate_dimensionless_params(espresso_params)
    extraction_percentage(sol, dimensionless_params, sim_params; cumul=cumul)
end

"""
    extraction_percentage(sol, dimensionless_params, sim_params; 
                            cumul=false)

**Arguments:**
- `sol::ODESolution`
- `dimensionless_params::DimensionlessParameters`
- `sim_params::SimulationParameters`
- `cumul::Bool=false`
"""
function extraction_percentage(sol:: ODESolution, dimensionless_params:: DimensionlessParameters, sim_params:: SimulationParameters; cumul::Bool=false)
    ϕs = dimensionless_params.ϕs
    q_z = dimensionless_params.q_z
    β = dimensionless_params.β
    Nx = sim_params.Nx

    if !cumul
        extract = trapz(sol.t, sol[Nx, :])
    else
        extract = cumul_integrate(sol.t, sol[Nx, :])
    end
    q_z*β*extract/ϕs
end

function extraction_yield(sol:: ODESolution, espresso_params:: EspressoParameters, sim_params:: SimulationParameters; cumul::Bool=false)
    dimensionless_params = generate_dimensionless_params(espresso_params)
    extraction_yield(sol, dimensionless_params, sim_params; cumul=cumul)
end

function extraction_yield(sol::ODESolution, dimensionless_params:: DimensionlessParameters, sim_params:: SimulationParameters; cumul::Bool=false)
    α = dimensionless_params.α
    α * extraction_percentage(sol, dimensionless_params, sim_params; cumul=cumul)
end