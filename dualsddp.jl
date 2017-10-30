

"""Init primal SDDP interface"""
function initprimal()
    model = build_model()
    params = getparams()
    sddpprimal = SDDPInterface(model, params, SDDP.IterLimit(MAX_ITER),
                               verbose_it=10)
    SDDP.init!(sddpprimal)
    return sddpprimal
end


"""Init dual SDDP interface"""
function initdual(sddp)
    modeldual = buildemptydual(sddp.spmodel.noises)
    sddpdual = SDDPInterface(modeldual, sddp.params,
                            SDDP.IterLimit(MAX_ITER),
                            verbose_it=0)
    initdual!(sddpdual)
    return sddpdual
end


"""Run SDDP primal alone."""
function runprimal!(sddpprimal)
    println("RUN PRIMAL SDDP")
    ubp = []
    stdp = []
    scen = SDDP.simulate_scenarios(sddpprimal.spmodel.noises, MCSIZE)

    tic()
    for iter in 1:MAXIT
        SDDP.iteration!(sddpprimal)

        if false #iter % UPPER_BOUND == 0
            cost = SDDP.simulate(sddpprimal, scen)[1]
            push!(ubp, mean(cost))
            push!(stdp, std(cost))
        end
        (iter % 10 == 0) && SDDP.reload!(sddpprimal)
    end
    texec = toq()
    println("Primal exec time: ", texec)
    SAVE && writecsv("lbprimal", sddpprimal.stats.lower_bounds)
    SAVE && writecsv("timeprimal", sddpprimal.stats.exectime)
end


"""Run SDDP dual alone"""
function rundual!(sddpdual, sddpprimal)
    ubd = []
    stdd = []
    srand(2713)
    scen = SDDP.simulate_scenarios(sddpdual.spmodel.noises, MCSIZE)

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
        (iter % 10 == 0) && displayit(iter, lb)
    end
    texec = toq()
    println("Dual exec time: ", texec)
    SAVE && writecsv("lbdual", lbdual)
    SAVE && writecsv("timedual", sddpdual.stats.exectime)
    return lbdual
end


function runjoint!(sddpprimal, sddpdual)
    lbdual = Float64[]
    timedual = Float64[]
    ubp = []
    stdp = []
    scen = SDDP.simulate_scenarios(sddpdual.spmodel.noises, MCSIZE)
    println("RUN DUAL SDDP")
    tic()

    for iter in 1:MAXIT
        td = SDDP.iteration!(sddpprimal, sddpdual)
        tic()
        lb = updateinitialstate!(sddpdual, X0)
        SDDP.fwdcuts(sddpdual)
        lb = updateinitialstate!(sddpdual, X0)
        tdual = toq() + td

        # save current iterations
        push!(lbdual, lb)
        push!(timedual, tdual)
        if iter % 50 == 0
            cost = SDDP.simulate(sddpprimal, scen)[1]
            push!(ubp, mean(cost))
            push!(stdp, std(cost))
        end

        (iter % 10 == 0) && displayit(iter, lb)
        (iter % 10 == 0) && SDDP.reload!(sddpprimal)
    end
    texec = toq()
    println("Total exec time: ", texec)

    return lbdual, timedual, ubp, stdp

end
