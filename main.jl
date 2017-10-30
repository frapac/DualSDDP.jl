# Run dual SDDP on dual of dam1

include("utils.jl")
include("config.jl")
include("dualutils.jl")
include("dualsddp.jl")
include("mpts2.jl")
include("innerapprox.jl")


# params
SAVE   = false
MAXIT  = 1000
NSIMU  = 1000
MCSIZE = 1000
PRIMAL = true
DUAL   = false

sddpprimal = initprimal()
sddpdual = initdual(sddpprimal)

lbdual, timedual, ubp, stdp = runjoint!(sddpprimal, sddpdual)
timedual -= sddpprimal.stats.exectime


### RESULTS
lbprimal = sddpprimal.stats.lower_bounds[end]
ubdual = lbdual[end]

println("#"^70)
println("Results --- $MAXIT iterations")
println("-------")
println("Primal LB:\t", lbprimal)
println("Dual UB:\t", ubdual)
println("Gap:\t", abs(lbprimal-ubdual)/lbprimal)


