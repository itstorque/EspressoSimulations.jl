Base.@kwdef mutable struct EspressoParameters

    R0:: Float64
    L:: Float64
    a1:: Float64
    a2:: Float64
    ρ_out:: Float64
    ρ_grounds:: Float64
    csat:: Float64
    cs0:: Float64
    ϕs:: Float64
    Mout:: Float64
    P:: Float64
    Ds_star:: Float64
    Deff_star:: Float64
    k:: Float64
    tshot:: Float64
    vol_fracs:: Float64

end

function Base.:*(a::Number, b::EspressoParameters)::EspressoParameters
    EspressoParameters([a * getproperty(b, fname) for fname in fieldnames(EspressoParameters)]...)
end

function Base.:*(a::EspressoParameters, b::EspressoParameters)::EspressoParameters
    EspressoParameters([getproperty(a, fname) * getproperty(b, fname) for fname in fieldnames(EspressoParameters)]...)
end

function Base.:+(a::EspressoParameters, b::EspressoParameters)::EspressoParameters
    EspressoParameters([getproperty(a, fname) + getproperty(b, fname) for fname in fieldnames(EspressoParameters)]...)
end

function EspressoParameters(espresso_params:: EspressoParameters; kwargs...)
    new_params = copy(espresso_params)
    for (i, j) in pairs(kwargs)
        setproperty!(new_params, i, j)
    end
    new_params
end

Base.copy(s::EspressoParameters) = EspressoParameters(
    R0=s.R0,
    L=s.L,
    a1=s.a1,
    a2=s.a2,
    ρ_out=s.ρ_out,
    ρ_grounds=s.ρ_grounds,
    csat=s.csat,
    cs0=s.cs0,
    ϕs=s.ϕs,
    Mout=s.Mout,
    P=s.P,
    Ds_star=s.Ds_star,
    Deff_star=s.Deff_star,
    k=s.k,
    tshot=s.tshot,
    vol_fracs=s.vol_fracs,
)

struct SimulationParameters
    Nx:: Int# = 40
    dx:: Float64# = 1 / (N - 1)
    x::LinRange# = LinRange(0, 1, N)

    T_end:: Float64# = 10

    function SimulationParameters(; Nx:: Int, T_end:: Float64)
        dx = 1/(Nx-1)
        x = LinRange(0, 1, Nx)
        return new(Nx, dx, x, T_end)
    end

end

Base.@kwdef struct DimensionlessParameters

    bet1:: Float64
    bet2:: Float64
    Ds1:: Float64
    Ds2:: Float64
    Q1:: Float64
    Q2:: Float64
    β:: Float64
    q_z:: Float64  # z-component of darcy flux
    Deff:: Float64
    K:: Float64
    α:: Float64
    ϕs:: Float64

end

function generate_dimensionless_params(espresso_params::EspressoParameters)
    
    R0 = espresso_params.R0
    L = espresso_params.L
    a1 = espresso_params.a1
    a2 = espresso_params.a2
    ρ_out = espresso_params.ρ_out
    ρ_grounds = espresso_params.ρ_grounds
    csat = espresso_params.csat
    cs0 = espresso_params.cs0
    ϕs = espresso_params.ϕs
    Mout = espresso_params.Mout
    P = espresso_params.P
    Ds_star = espresso_params.Ds_star
    Deff_star = espresso_params.Deff_star
    k = espresso_params.k
    tshot = espresso_params.tshot

    # The vol. fracs of fines and boulders (of the grounds)
    f1 = espresso_params.vol_fracs
    f2 = 1 - f1

    # Calculate vol. fracs. of fines and boulders (of baskets)
    ϕs1 = ϕs * f1
    ϕs2 = ϕs * f2

    # Brunauer-Emmett-Teller surface area of particles with radius a1 and a2
    # in 1/m units
    bet1_star = 3 * ϕs1 / a1
    bet2_star = 3 * ϕs2 / a2
    bet0_star = bet1_star / 2 + bet2_star / 2

    # Convert to dimensionless parameters
    bet1 = bet1_star / bet0_star
    bet2 = bet2_star / bet0_star
    Ds1 = Ds_star * tshot / a1^2
    Ds2 = Ds_star * tshot / a2^2
    Q1 = 1 / a1 / bet0_star
    Q2 = 1 / a2 / bet0_star
    β = csat / cs0
    q_z = Mout / π / R0^2 / ρ_out / L
    α = cs0 / ρ_grounds
    K = k * cs0^2 * tshot * bet0_star
    Deff = Deff_star * tshot / L^2

    return DimensionlessParameters(
        bet1=bet1,
        bet2=bet2,
        Ds1=Ds1,
        Ds2=Ds2,
        Q1=Q1,
        Q2=Q2,
        β=β,
        q_z=q_z,
        Deff=Deff,
        K=K,
        α=α,
        ϕs=ϕs
    )
end

function model_params_to_vec(espresso_params:: EspressoParameters, sim_params::SimulationParameters)
    dimensionless_params = generate_dimensionless_params(espresso_params)
    model_params_to_vec(dimensionless_params, sim_params)
end

function model_params_to_vec(dimensionless_params:: DimensionlessParameters, sim_params::SimulationParameters)

    Nx, dx, x = sim_params.Nx, sim_params.dx, sim_params.x

    Deff = dimensionless_params.Deff
    Ds1 = dimensionless_params.Ds1
    Ds2 = dimensionless_params.Ds2
    bet1 = dimensionless_params.bet1
    bet2 = dimensionless_params.bet2
    K = dimensionless_params.K
    Q1 = dimensionless_params.Q1
    Q2 = dimensionless_params.Q2
    β = dimensionless_params.β
    ϕs = dimensionless_params.ϕs
    q_z = dimensionless_params.q_z

    (Nx, dx, x, Deff, Ds1, Ds2, bet1, bet2, K, Q1, Q2, β, ϕs, q_z)

end