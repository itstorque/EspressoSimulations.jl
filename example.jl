using EspressoSimulations
using DifferentialEquations
using Plots

sim_params = SimulationParameters(Nx=40, T_end=10.)

espresso_params = EspressoParameters(
    R0 = 29.2e-3, # puck radius
    L = 18.7e-3, # Cylinder height
    a1 = 80e-6, # fine grain radii [m]
    a2 = 300e-6, # coarse grain radii [m]
    ρ_out = 997, # density of water at 90 deg C [kg/m3]
    ρ_grounds = 330, # density of roasted coffee bulk [kg/m3]
    csat = 2.124e2, # saturation concentration of water outside grain. temperature dependent [kg/m3]
    cs0 = 1.18e2, # concentration of solubles in the grains initially [kg/m3], assuming preinfusion, see https:://researchrepository.ul.ie/articles/thesis/Heat_and_mass_transfer_in_dispersed_two-phase_flows/19814533/1?file=35260576 for pre-infusion details
    ϕs = 0.8272, # volume fraction of grounds in a packed bed [unitless] (http:://refhub.elsevier.com/S2590-2385(19)30410-2/sref30)
    Mout = 0.04, # mass of the beverage [kg]
    P = 5, # pressure [bar]
    Ds_star = 0.625e-9, # 
    Deff_star = 1e-6,
    k = 1e-9, # reaction rate constant
    tshot = 33.9,

    vol_fracs = 0.1624,
);

ep = copy(espresso_params)

sweep_ϕs = [0.5, 0.7, 0.8]
sols = Vector{ODESolution}(undef, length(sweep_ϕs))

for (iϕs, ϕs) in enumerate(sweep_ϕs)
    ep.ϕs = ϕs
    prob = to_odeproblem(ep, sim_params)
    sols[iϕs] = solve(prob, Rodas5(autodiff=false))
end

plot(legendtitle="\$\\phi s\$")
ylabel!("EY Percentage")
xlabel!("time")

for (ϕs, sol) in zip(sweep_ϕs, sols)
    plot!(sol.t, extraction_percentage(sol, espresso_params, sim_params; cumul=true); label=ϕs)
end

plot!()