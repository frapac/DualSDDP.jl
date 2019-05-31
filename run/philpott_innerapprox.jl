srand(1111)

include("src/dualSDDP.jl")

# define a production transport problem MPTS
mpts = MPTS([:FRA, :GER, :ESP, :UK, :PT, :ITA, :SUI, :BEL], 12, 10)

# Init SDDP interface
sddpprimal = initprimal(mpts)

ubp, stdp,  trajs = runprimal!(sddpprimal, maxiterations=1000)

include("src/dp.jl")
println(solvedp!(sddpprimal.spmodel, trajs)[1,1])
