using LinearAlgebra
using SparseArrays

function build_mass(sim_params::SimulationParameters)

    N, dr, r = sim_params.Nx, sim_params.dx, sim_params.x

    @assert r[1] == 0
    @assert r[N] == 1
    @assert length(r) == N

    Md1 = 4π * (r.^2 * dr .+ dr^3/12)

    Md1[1] = (4/3) * π * (dr/2)^3
    Md1[N] = 4/3 * π - 4/3 * π * (1 - dr/2)^3

    M1 = Diagonal(Md1)

    M2 = Tridiagonal(
        1/8*ones(N-1),
        6/8*ones(N),
        1/8*ones(N-1)
    )

    M2[:, 1] *= 2
    M2[:, N] *= 2
    M2[1, 1] = 3/4
    M2[N, N] = 3/4
    
    M_sub = M2 * M1
    
    # Create the full mass matrix
    M_size = (2*N+1)*N
    M = spzeros(Float64, M_size, M_size)

    for i in 1:2*N
        M[i*N+1:i*N+N, i*N+1:i*N+N] .= M_sub
    end
    
    # Set algebraic boundary conditions
    M[1, 1] = 0
    M[N, N] = 0
    
    M
end

function init_problem_vector(sim_params::SimulationParameters)

    cl0 = zeros(sim_params.Nx) # initial concentration of sollubles in liquid
    u10 = ones(sim_params.Nx)
    u20 = ones(sim_params.Nx)

    # Concatenate the initial data
    u0 = vcat(cl0, u10)
    for _ in 1:sim_params.Nx-1
        u0 = vcat(u0, u10)
    end
    
    for _ in 1:sim_params.Nx
        u0 = vcat(u0, u20)
    end
    
    u0

end

function RHS(du, u, params, t)

    (N, dx, x, Deff, Ds1, Ds2, bet1, bet2, K, Q1, Q2, β, ϕs, q_z) = params

    cl = @view u[1:N]
    u1 = @view u[N+1:N*(N+1)]
    u2 = @view u[N*(N+1)+1:end]

    dcl = @view du[1:N]
    du1 = @view du[N+1:N*(N+1)]
    du2 = @view du[N*(N+1)+1:end]

    ru1 = reshape(u1, (N, N))
    rdu1 = reshape(du1, (N, N))
    ru2 = reshape(u2, (N, N))
    rdu2 = reshape(du2, (N, N))

    # Get surface concentrations
    csurf1 = ru1[N, :]
    csurf2 = ru2[N, :]

    # Compute the reaction rate
    G1 = K * (1 .- cl) .* csurf1 .* (csurf1 - β .* cl)
    G2 = K * (1 .- cl) .* csurf2 .* (csurf2 - β .* cl)

    # Eqn for c_l (concentration in the liquid)
    gamma = 1 / (1 - ϕs)

    dcl[1] = -Deff / dx * (-3cl[1]/2 + 2cl[2] - cl[3]/2) + q_z * cl[1]
    
    dcl[2:N-1] .= -q_z * gamma * (cl[(2:N-1).+1] .- cl[(2:N-1).-1])/2dx .+
                    gamma * Deff * (cl[(2:N-1).+1] .- 2cl[(2:N-1)] .+ cl[(2:N-1).-1])/dx^2 .+
                    gamma * bet1 * G1[(2:N-1)] .+
                    gamma * bet2 * G2[(2:N-1)]

    dcl[N] = cl[N-2]/2 - 2cl[N-1] + 3cl[N]/2

    # Eqns for cs1
    rad_window = 1:N-1

    t1 = 4π * (x[rad_window.+1]/2 + x[rad_window]/2).^2 .* Ds1 .* (ru1[rad_window.+1, :] .- ru1[rad_window, :])/dx

    rdu1[1:N-1, :] = t1
    rdu1[N, :] = -4π * β * Q1 * G1 # set boundary condition

    rdu1[2:N, :] .-= t1
    

    # Eqns for cs2
    t2 = 4π * (x[rad_window.+1]/2 .+ x[rad_window]/2).^2 .* Ds2 .* (ru2[rad_window.+1, :] .- ru2[rad_window, :])/dx

    rdu2[1:N-1, :] = t2
    rdu2[N, :] = -4π * β * Q2 * G2 # set boundary condition

    rdu2[2:N, :] .-= t2


end

function to_odeproblem(problem_params:: Union{EspressoParameters, DimensionlessParameters}, sim_params:: SimulationParameters)
    
    M = build_mass(sim_params)
    u0 = init_problem_vector(sim_params)

    f = ODEFunction(RHS, mass_matrix=M)
    prob = ODEProblem(f, u0, (0.0, sim_params.T_end), model_params_to_vec(problem_params, sim_params));
end