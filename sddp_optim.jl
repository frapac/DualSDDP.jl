
################################################################################
################################################################################
# DAMS VALLEY MANAGEMENT
# -----------
# This file specify the parameters of the problem
################################################################################
################################################################################

using JuMP, StochDynamicProgramming, Clp, Gurobi
const SDDP = StochDynamicProgramming

# solver
using Gurobi
SOLVER = Gurobi.GurobiSolver(OutputFlag=false, Threads=1)

# problem
include("config.jl")
include("dualutils.jl")
include("damsvalley.jl")

srand(1111)

# Read costs along time:
const α = EFFICIENCY*COST_CONV
const COST = collect(Float64, readdlm("$DATA_SOURCE/generated_prices.txt") * α)

# Read intakes:
INTAKES = [readdlm("$DATA_SOURCE/generated_intakes_$i.txt") for i in 0:N_DAMS-1]


# Build noiselaws with input data:
function build_noiselaws()
    laws = Vector{SDDP.NoiseLaw}(TF-1)
    n_scenarios = size(INTAKES[1], 1)
    # Uniform probabilities:
    proba = 1./n_scenarios*ones(n_scenarios)
    for t in 1:TF-1
        alea = hcat([INTAKES[i][:, t] for i in 1:N_DAMS]...)'
        laws[t] = SDDP.NoiseLaw(alea, proba)
    end
    return laws
end


# Build dynamic model and SDDP solver:
function build_model()
    srand(1111)
    laws = build_noiselaws();
    i = ones(Float64, N_DAMS)
    c = vcat(-i, 0*i)

    A, B, D, E, _ = getmatrix(N_DAMS, 1)
    dynamic(t, x, u, w) = x + B*u + w
    cost_t(t, x, u, w) = COST[t] * dot(c, u)

    # Build bounds:
    x_bounds = [(STOCK_MIN, STOCK_MAX) for i in 1:N_DAMS];
    u_bounds = vcat([(CONTROL_MIN, CONTROL_MAX) for i in 1:N_DAMS], [(0., Inf) for i in 1:N_DAMS]);

    model = SDDP.LinearSPModel(TF,       # number of timestep
                               u_bounds, # control bounds
                               X0,       # initial state
                               cost_t,   # cost function
                               dynamic,  # dynamic function
                               laws,
                              )
    SDDP.set_state_bounds(model, x_bounds)

    params = SDDP.SDDPparameters(SOLVER, passnumber=FORWARD_PASS,
                                 gap=EPSILON, max_iterations=MAX_ITER,
                                 confidence=CONFIDENCE_LEVEL,
                                 montecarlo_in_iter=MONTE_CARLO_SIZE,
                                 montecarlo_final=FINAL_MONTE_CARLO_SIZE,
                                 prune_cuts=PRUNE_CUTS,
                                 pruning_algo="exact+",
                                 rho0=0, alpha=.95)
    params.compute_ub = UPPER_BOUND
    return model, params
end

