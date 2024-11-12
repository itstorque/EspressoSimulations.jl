using LinearAlgebra
using SparseArrays

"""
    slice_volumes(sim_params)

    Generate a list of slice volumes for a spherical diffusion problem, 
    where \$r\$ is the center coordinate and \$dr\$ is the length of the change in 
    radius. The change in radius is centered around \$r\$. The first and last slice
    have a shell of radius \$dr/2\$.

    Every shell has a volume of \$4π \\cdot ((r + dr)^2 - (r - dr)^2)\$.
"""
function slice_volumes(sim_params::SimulationParameters)
    N, dr, r = sim_params.Nx, sim_params.dx, sim_params.x

    # vol of slice i
    Vi = 4π * ((r .+ dr).^2 .- (r .- dr).^2)

    # vol of first slice
    Vi[1] = (4/3) * π * (dr/2)^3
    # vol of last slice
    Vi[N] = 4/3 * π - 4/3 * π * (1 - dr/2)^3

    Vi

end

"""
    mass_matrix(N, Vi)

    Generate a tile for a block mass matrix.
        
    This differs by a factor of (1-ϕs) from the publication's supplementary [1]
"""
function mass_matrix(N::Integer, Vi::Vector{Float64})
    # Equations 20 and 23 from https://dspace.mit.edu/bitstream/handle/1721.1/91234/Bazant_Efficient%20conservative.pdf?sequence=1&isAllowed=y
    # the volumes here are equidistant, need to generalize this later
    M1 = Tridiagonal(
        1/8*ones(N-1),
        6/8*ones(N),
        1/8*ones(N-1)
    )

    M1[:, 1] *= 2
    M1[:, N] *= 2
    M1[1, 1] = 3/4
    M1[N, N] = 3/4
    
    M2 = Diagonal(Vi)
    Mc = M1 * M2

    Mc
end

"""
    mass_matrix(sim_params)

    Build the mass matrix for our simulation from a set of 
    `SimulationParameters`. This will generate a \$2N^2+N\$
    block tridiagonal matrix. Block diagonal matrix using sub-blocks
    from `mass_matrix(N, Vi)`.

    This differs by a factor of (1-ϕs) from the publication's supplementary [1]

    **Arguments**
    - `sim_params::SimulationParameters`: simulation parameters to construct a mass matrix for
"""
function mass_matrix(sim_params::SimulationParameters)

    N = sim_params.Nx

    Vi = slice_volumes(sim_params)

    M_sub = mass_matrix(N, Vi)
    
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

"""
    init_problem_vector(sim_params)

Initial conditions for a simulation with the parameters `sim_params`. Initializes an
\$2N^2+N\$ vector where:
- the first \$N\$ entries represent the initial concentration of sollubles in the 
liquid \$c_l(t=0) = 0\$.
- The next \$N^2\$ entries encode a \$N\\times N\$ matrix that \$u_1(t=0) = 1\$
- The last \$N^2\$ entries encode a \$N\\times N\$ matrix that \$u_2(t=0) = 1\$

**Arguments**
- `sim_params::SimulationParameters`: simulation parameters to build initial conditions for.
"""
function init_problem_vector(sim_params::SimulationParameters)

    cl0 = zeros(sim_params.Nx) # initial concentration of sollubles in liquid
    u10 = ones(sim_params.Nx^2)
    u20 = ones(sim_params.Nx^2)

    # Concatenate the initial data
    vcat(cl0, u10, u20)

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
    # see (18) in paper
    gamma = 1 / (1 - ϕs)
    
    # equation 18
    dcl[2:N-1] .= -q_z * gamma * (cl[(2:N-1).+1] .- cl[(2:N-1).-1])/2dx .+
                    gamma * Deff * (cl[(2:N-1).+1] .- 2cl[(2:N-1)] .+ cl[(2:N-1).-1])/dx^2 .+
                    gamma * bet1 * G1[(2:N-1)] .+
                    gamma * bet2 * G2[(2:N-1)]

    # equations 18 + (19 & 20)
    dcl[1] = -Deff / dx * (-3cl[1]/2 + 2cl[2] - cl[3]/2) + q_z * cl[1]
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
    
    M = mass_matrix(sim_params)
    u0 = init_problem_vector(sim_params)

    f = ODEFunction(RHS, mass_matrix=M)
    prob = ODEProblem(f, u0, (0.0, sim_params.T_end), model_params_to_vec(problem_params, sim_params));
end