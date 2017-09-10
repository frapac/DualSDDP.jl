
################################################################################
################################################################################
# DAMS VALLEY MANAGEMENT
# -----------
# stochdec.jl
# -----------
# This file specify the parameters of the problem
################################################################################
################################################################################

import JuMP, StochDynamicProgramming, Clp, Gurobi
const SDDP = StochDynamicProgramming

include("config.jl")
include("problem.jl")
include("dual.jl")
include("modeler.jl")

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

    # Build bounds:
    x_bounds = X_BOUNDS
    u_bounds = U_BOUNDS
    # Define target for stock:
    x0 = STOCK_TARGET

    model = SDDP.LinearSPModel(TF,       # number of timestep
                               u_bounds, # control bounds
                               x0,       # initial state
                               cost_t,   # cost function
                               dynamic,  # dynamic function
                               laws,
                               #= Vfinal=final_cost_dams, =#
                               #= eqconstr=constdual =#
                              )
    SDDP.set_state_bounds(model, x_bounds)

    solver = Gurobi.GurobiSolver(OutputFlag=false, Threads=1)
    params = SDDP.SDDPparameters(solver, passnumber=FORWARD_PASS,
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

