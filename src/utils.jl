################################################################################
################################################################################
# SDDP DUAL
# -----------
# This file specify the parameters of the problem
################################################################################
################################################################################

 SOLVER = Gurobi.GurobiSolver(OutputFlag=false, Threads=1) 
#= SOLVER = Xpress.XpressSolver(PRESOLVE=0, OUTPUTLOG=2) =#
#SOLVER = Clp.ClpSolver(LogLevel=0, PresolveType=0)
#= Xpress.setparam!(SOLVER, Xpress.XPRS_PRESOLVE,0) =#



"Get SDDP parameters for StochDynamicProgramming."
function getparams(;maxit=50, nbsimu=1000, Δub=100, Δprune=100)
    params = SDDP.SDDPparameters(SOLVER, passnumber=1,
                                 gap=0.001, max_iterations=maxit,
                                 confidence=.975,
                                 montecarlo_in_iter=nbsimu,
                                 montecarlo_final=nbsimu,
                                 prune_cuts=Δprune,
                                 pruning_algo="exact+",
                                 rho0=0, alpha=.95)
    params.compute_ub = Δub
    return params
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
