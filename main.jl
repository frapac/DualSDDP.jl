################################################################################
# CERMICS, ENPC
# SDDP dual
# main code
################################################################################

srand(1111)

include("src/dualSDDP.jl")

# params
SAVE   = false
# 1: Primal SDDP   2: Dual SDDP    3: Mix primal / dual
# 4: computing Philpott's upper bound
ALGO = 4

# define a production transport problem MPTS
mpts = MPTS([:FRA, :GER, :ESP, :UK, :PT, :ITA, :SUI, :BEL], 36, 10)

# Init SDDP interface
sddpprimal = initprimal(mpts)

if ALGO == 1
    ubp, stdp = runprimal!(sddpprimal)
elseif ALGO == 2
    sddpdual = initdual(mpts, sddpprimal)
    lbdual, timedual = rundual!(sddpdual, sddpprimal)
elseif ALGO == 3
    sddpdual = initdual(mpts, sddpprimal)
    lbdual, timedual, ubp, stdp = runjoint!(sddpprimal, sddpdual, maxiterations=100)
    # recalibrate dual time
    timedual -= sddpprimal.stats.exectime
elseif ALGO == 4
    include("src/dp.jl")
    ubp, stdp,  trajs = runprimal!(sddpprimal, maxiterations=1000)
    ubphilpott = []
    for n in [50,100,200,300,400,500,600,700,800,900,1000]
        println(n)
        ub_p =  solvedp!(sddpprimal.spmodel, trajs[1:n])[1,1]
        println(ub_p)
        push!(ubphilpott , ub_p )
    end
end


### RESULTS
#lbprimal = sddpprimal.stats.lower_bounds[end]

# println("#"^70)
# println("Results")
# println("-------")
# println("Primal LB:\t", lbprimal)
# if ALGO >= 2
#     ubdual = lbdual[end]
#     println("Dual UB:\t", ubdual)
#     println("Gap:\t", abs(lbprimal-ubdual)/lbprimal)
# end
