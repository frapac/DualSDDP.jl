################################################################################
# CERMICS, ENPC
# SDDP dual
################################################################################
# optimization packages
using JuMP, StochDynamicProgramming, Clp
# use Gurobi as LP solver
using Gurobi
# import Clustering for kmeans quantization of noises
using Clustering

const SDDP = StochDynamicProgramming

# utililities
include("utils.jl")
include("dualutils.jl")
include("dualsddp.jl")
include("data.jl")
include("mpts2.jl")
include("innerapprox.jl")


# params
SAVE   = false
MAXIT  = 1000
NSIMU  = 1000
MCSIZE = 1000
# 1: Primal SDDP   2: Dual SDDP    3: Mix primal / dual
ALGO = 1
PRIMAL = true
DUAL   = false
ΔMC = 100

# Init SDDP interface
sddpprimal = initprimal()

if ALGO == 1
    ubp, stdp = runprimal!(sddpprimal)
elseif ALGO == 2
    sddpdual = initdual(sddpprimal)
    lbdual, timedual = rundual!(sddpdual, sddpprimal)
elseif ALGO == 3
    sddpdual = initdual(sddpprimal)
    lbdual, timedual, ubp, stdp = runjoint!(sddpprimal, sddpdual)
    # recalibrate dual time
    timedual -= sddpprimal.stats.exectime
end


### RESULTS
lbprimal = sddpprimal.stats.lower_bounds[end]
ubdual = lbdual[end]

println("#"^70)
println("Results --- $MAXIT iterations")
println("-------")
println("Primal LB:\t", lbprimal)
println("Dual UB:\t", ubdual)
println("Gap:\t", abs(lbprimal-ubdual)/lbprimal)
