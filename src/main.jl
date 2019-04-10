################################################################################
# CERMICS, ENPC
# SDDP dual
# main code
################################################################################

srand(1111)

include("dualSDDP.jl")

# params
SAVE   = false
# 1: Primal SDDP   2: Dual SDDP    3: Mix primal / dual
ALGO = 1

# define a production transport problem MPTS
mpts = MPTS([:FRA, :GER], 5, 10)

# Init SDDP interface
sddpprimal = initprimal(mpts)

if ALGO == 1
    ubp, stdp, t = runprimal!(sddpprimal)
elseif ALGO == 2
    sddpdual = initdual(mpts, sddpprimal)
    lbdual, timedual = rundual!(sddpdual, sddpprimal)
elseif ALGO == 3
    sddpdual = initdual(mpts, sddpprimal)
    lbdual, timedual, ubp, stdp = runjoint!(sddpprimal, sddpdual)
    # recalibrate dual time
    timedual -= sddpprimal.stats.exectime
end


### RESULTS
lbprimal = sddpprimal.stats.lower_bounds[end]

println("#"^70)
println("Results")
println("-------")
println("Primal LB:\t", lbprimal)
if ALGO >= 2
    ubdual = lbdual[end]
    println("Dual UB:\t", ubdual)
    println("Gap:\t", abs(lbprimal-ubdual)/lbprimal)
end
