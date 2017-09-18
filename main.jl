# Run dual SDDP on dual of dam1

include("sddp_optim.jl")
include("innerapprox.jl")

# params
MAXIT = 50
NSIMU = 1000
PRIMAL = true
DUAL = true

# OUTER APPROXIMATION
OA = true
# INNER APPROXIMATION
IA = false
# JOINT APPROXIMATION
JA = false

### Build primal problem
model, params = build_model()
sddpprimal = SDDPInterface(model, params,
                    SDDP.IterLimit(MAX_ITER),
                    verbose_it=10)
SDDP.init!(sddpprimal)

### Build dual problem
modeldual = buildual_dam1(build_noiselaws())
sddpdual = SDDPInterface(modeldual, params,
                         SDDP.IterLimit(MAX_ITER),
                         verbose_it=0)
initdual!(sddpdual)


### SDDP DUAL ####
if DUAL
    # Run 1 combined iteration to init cuts in sddpdual
    for iter in 1:1
        SDDP.iteration!(sddpprimal, sddpdual)
    end
    p0 = SDDP.get_subgradient(sddpprimal.bellmanfunctions[1], X0)


    ### RUN iterations in dual
    lbdual = Float64[]
    println("RUN DUAL SDDP")
    lb = updateinitialstate!(sddpdual, X0)
    tic()
    for iter in 1:MAXIT
        # Update initial costate
        lb = updateinitialstate!(sddpdual, X0)

        # Run forward an backward pass
        SDDP.iteration!(sddpdual)

        # save current iterations
        push!(lbdual, lb)
        (iter % 10 == 0) && display(iter, lb)
    end
    texec = toq()
    println("Dual exec time: ", texec)
    println(sddpdual.spmodel.initialState)
end


### RUN iterations in primal
if PRIMAL
    println("RUN PRIMAL SDDP")
    tic()
    for iter in 1:MAXIT
        SDDP.iteration!(sddpprimal)
    end
    texec = toq()
    println("Primal exec time: ", texec)
end


### MONTE CARLO ESTIMATION
if OA
    println("RUN OUTER APPROX")
    srand(2713)
    c = @time SDDP.simulate(sddpprimal, NSIMU)[1]
end


### INNER APPROX
if IA
    println("RUN INNER APPROX")
    # Build a new SDDP interface
    sddp= SDDPInterface(model, params,
                        SDDP.IterLimit(MAX_ITER),
                        verbose_it=10)
    # Replace model with inner approx
    init_innermodeler!(sddp, sddpdual.bellmanfunctions)
    #= println(sddp.solverinterface[2]) =#
    srand(2713)
    ci = @time SDDP.simulate(sddp, NSIMU)[1]
end


### JOINT APPROX
if JA
    println("RUN JOINT APPROX")
    # Build a new SDDP interface
    jointcost = Float64[]
    for w in 0:.1:1
        sddp= SDDPInterface(model, params,
                            SDDP.IterLimit(MAX_ITER),
                            verbose_it=10)
        # Replace model with inner approx
        init_jointmodeler!(sddp, sddpdual.bellmanfunctions, sddpprimal.bellmanfunctions, w)
        srand(2713)
        cj = @time SDDP.simulate(sddp, NSIMU)[1]
        push!(jointcost, mean(cj))
    end
end


### RESULTS
lbprimal = sddpprimal.stats.lower_bounds[end]
ubdual = lbdual[end]

println("#"^70)
println("Results")
println("-------")
println("Primal LB:\t", lbprimal)
println("Dual UB:\t", -ubdual)
println("Gap:\t", (lbprimal+ubdual)/lbprimal)
OA && println("Monte Carlo (OA):\t", mean(c))
IA && println("Monte Carlo (IA):\t", mean(ci))
JA && println("Monte Carlo (JA):\t", jointcost)

