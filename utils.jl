
################################################################################
################################################################################
# SDDP DUAL
# -----------
# This file specify the parameters of the problem
################################################################################
################################################################################

SOLVER = Gurobi.GurobiSolver(OutputFlag=false, Threads=1)

# problem
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

"""Build incidence matrix and return flows' bounds."""
function buildincidence{T}(connexion::Matrix{T})
    nnodes = size(connexion, 1)
    narcs = floor(Int, sum(connexion .> 0.)/2)

    # incidence matrix
    A = zeros(Int, nnodes, narcs)

    ic = 0
    bounds = Float64[]
    for ix in 1:(nnodes-1), iy in (ix+1):nnodes
        if connexion[ix, iy] > 0
            ic += 1
            A[ix, ic] =  1
            A[iy, ic] = -1
            push!(bounds, connexion[ix, iy])
        end
    end
    return A, bounds
end
