
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
include("mpts.jl")

srand(1111)


function getparams()
    params = SDDP.SDDPparameters(SOLVER, passnumber=FORWARD_PASS,
                                 gap=EPSILON, max_iterations=MAX_ITER,
                                 confidence=CONFIDENCE_LEVEL,
                                 montecarlo_in_iter=MONTE_CARLO_SIZE,
                                 montecarlo_final=FINAL_MONTE_CARLO_SIZE,
                                 prune_cuts=PRUNE_CUTS,
                                 pruning_algo="exact+",
                                 rho0=0, alpha=.95)
    params.compute_ub = UPPER_BOUND
    params
end
